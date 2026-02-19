#' ============================================================================
#' Robust Signal Resampling Function 
#' ============================================================================
#' 
#' A generic function for resampling signals to a target length using
#' polyphase filtering (via gsignal package).
#' 
#' @author Philippe Terrier / ORD 2025 Project (RACING)
#' @date February 2026
#' @version 2.1
#' @license MIT
#' 
#' This function handles the "close length issue" when the resampling 
#' ratio is close to 1, using two-stage resampling.
#' 
#' @references
#' MATLAB resample documentation:
#' https://www.mathworks.com/help/signal/ref/resample.html
#' ============================================================================

# =============================================================================
# Package Availability Check
# =============================================================================

# Check for gsignal (preferred)
HAS_GSIGNAL_PKG <- requireNamespace("gsignal", quietly = TRUE)

# Check for signal (fallback)
HAS_SIGNAL_PKG <- requireNamespace("signal", quietly = TRUE)

# Validate that at least one package is available
if (!HAS_GSIGNAL_PKG && !HAS_SIGNAL_PKG) {
  warning("Neither 'gsignal' nor 'signal' package is available.\n",
          "Install gsignal for MATLAB-compatible resampling:\n",
          "  install.packages('gsignal')\n",
          "Spline interpolation will be used as fallback (not recommended).")
}

if (!HAS_GSIGNAL_PKG && HAS_SIGNAL_PKG) {
  message("Note: Using 'signal' package. For better MATLAB compatibility,\n",
          "      install 'gsignal': install.packages('gsignal')")
}


# =============================================================================
# Internal Resampling Functions
# =============================================================================

#' Internal: Polyphase resampling using gsignal package #' 
#' gsignal::resample uses the same polyphase FIR algorithm as MATLAB:
#' 1. Design anti-aliasing FIR filter with Kaiser window
#' 2. Apply via upfirdn (polyphase decomposition)
#' 
#' @keywords internal
.resample_polyphase <- function(x, p, q) {
  if (!HAS_GSIGNAL_PKG) {
    stop("gsignal package not available")
  }
  return(gsignal::resample(x, p, q))
}


#' Internal: Bandlimited resampling using signal package (legacy fallback)
#' 
#' signal::resample uses sinc-based interpolation, which differs from MATLAB.
#' Use only as fallback when gsignal is not available.
#' 
#' @keywords internal
.resample_bandlimited <- function(x, p, q) {
  if (!HAS_SIGNAL_PKG) {
    stop("signal package not available")
  }
  return(signal::resample(x, p, q))
}


#' Internal: Spline-based resampling (last resort fallback)
#' 
#' Uses cubic spline interpolation. NOT recommended for dynamical analysis
#' as it may not preserve frequency content properly.
#' 
#' @keywords internal
.resample_spline <- function(x, target_length) {
  n <- length(x)
  x_orig <- seq(0, 1, length.out = n)
  x_new <- seq(0, 1, length.out = target_length)
  return(spline(x_orig, x, xout = x_new, method = "natural")$y)
}


# =============================================================================
# Utility Functions
# =============================================================================

#' Calculate Greatest Common Divisor (GCD)
#' 
#' Uses Euclidean algorithm to compute the GCD of two integers.
#' 
#' @param a First integer
#' @param b Second integer
#' @return GCD of a and b
#' 
#' @examples
#' gcd(18750, 16000)  # Returns 250
#' gcd(18750, 18100)  # Returns 50 (problematic - small GCD)
#' 
#' @export

gcd <- function(a, b) {
  a <- abs(as.integer(a))
  b <- abs(as.integer(b))
  while (b != 0) {
    temp <- b
    b <- a %% b
    a <- temp
  }
  return(a)
}


#' Check if Resampling Ratio is Problematic
#' 
#' Analyzes the resampling ratio to determine if it will cause
#' numerical instability in the polyphase filter.
#' 
#' @param original_length Integer. Length of the original signal.
#' @param target_length Integer. Desired length after resampling.
#' @param threshold Integer. Maximum allowed filter phases (default: 100).
#' 
#' @return A list containing:
#' \describe{
#'   \item{ratio}{The resampling ratio (target/original)}
#'   \item{gcd}{Greatest common divisor of lengths}
#'   \item{p_reduced}{Reduced numerator (target / gcd)}
#'   \item{q_reduced}{Reduced denominator (original / gcd)}
#'   \item{filter_phases}{Number of polyphase filter phases required}
#'   \item{is_problematic}{Logical. TRUE if ratio is problematic.}
#' }
#' 
#' @examples
#' # Safe ratio (large GCD)
#' check_resampling_ratio(16000, 18750)
#' # $filter_phases = 75, $is_problematic = FALSE
#' 
#' # Problematic ratio (small GCD, close to 1)
#' check_resampling_ratio(18100, 18750)
#' # $filter_phases = 375, $is_problematic = TRUE
#' 
#' @export
check_resampling_ratio <- function(original_length, target_length, threshold = 100) {
  

  # Input validation
  stopifnot(
    is.numeric(original_length), length(original_length) == 1, original_length > 0,
    is.numeric(target_length), length(target_length) == 1, target_length > 0,
    is.numeric(threshold), length(threshold) == 1, threshold > 0
  )
  
  original_length <- as.integer(original_length)
  target_length <- as.integer(target_length)
  
  # Calculate GCD
  g <- gcd(original_length, target_length)
  
  # Reduce to coprime integers
  p_reduced <- target_length / g
  q_reduced <- original_length / g
  
  # Calculate filter complexity (number of phases)
  # This is max(p, q) after reduction to coprime form
  filter_phases <- max(p_reduced, q_reduced)
  
  # Determine if problematic
  is_problematic <- filter_phases > threshold
  
  return(list(
    ratio = target_length / original_length,
    gcd = g,
    p_reduced = p_reduced,
    q_reduced = q_reduced,
    filter_phases = filter_phases,
    is_problematic = is_problematic
  ))
}


# =============================================================================
# Main Resampling Function
# =============================================================================

#' Resample Signal with Automatic Problematic Ratio Handling
#' 
#' Resamples a signal to a target length using polyphase filtering.
#' 
#' Automatically detects and handles problematic ratios where the
#' target length is close to the original length.
#' 
#' @param signal Numeric vector. The input signal to resample.
#' @param target_length Integer. The desired output length.
#' @param method Character. Resampling method:
#' \describe{
#'   \item{"auto"}{(default) Automatically detect and handle problematic ratios}
#'   \item{"direct"}{Always use direct resampling (may fail for problematic ratios)}
#'   \item{"two_stage"}{Always use two-stage resampling}
#' }
#' @param complexity_threshold Integer. Maximum filter phases before considering
#'   ratio problematic (default: 100).
#' @param intermediate_factor Numeric. Factor to multiply original length
#'   for intermediate step in two-stage resampling (default: 1.5 = 3/2).
#' @param debug Logical. Print diagnostic messages (default: FALSE).
#' @param use_fallback Logical. If TRUE and gsignal is unavailable,
#'   fall back to signal package, then spline (default: TRUE).
#' 
#' @return Numeric vector of length \code{target_length}. The resampled signal.
#'   Includes attributes:
#' \describe{
#'   \item{method_used}{"direct", "two_stage", or "none"}
#'   \item{original_length}{Length of input signal}
#'   \item{filter_phases}{Number of polyphase filter phases}
#'   \item{was_problematic}{Whether the ratio was detected as problematic}
#'   \item{backend}{"gsignal" (polyphase), "signal" (sinc), or "spline"}
#' }
#' 
#' @details
#' ## Algorithm (gsignal - MATLAB like)
#' 
#' The function uses polyphase filtering via \code{gsignal::resample()},
#' which implements the same algorithm as MATLAB's \code{resample()}:
#' 
#' 1. Design FIR anti-aliasing filter (Kaiser window)
#' 2. Apply via polyphase decomposition (upfirdn)
#' 3. Compensate for filter delay
#' 
#' This preserves the signal's frequency content up to the Nyquist frequency,
#' which is essential for preserving dynamical properties needed for:
#' - Lyapunov exponent estimation (LDS, ACI)
#' - Phase space reconstruction
#' - Attractor analysis
#' 
#' ## Problematic Ratio Handling
#' 
#' When the resampling ratio is close to 1, the polyphase filter requires
#' many phases, causing numerical instability. This is detected when
#' \code{max(p, q) / gcd(p, q) > complexity_threshold}.
#' 
#' The two-stage solution:
#' 1. Upsample by \code{intermediate_factor} (default 1.5)
#' 2. Resample to target length
#' 
#' This changes the GCD structure to reduce filter complexity.
#' 
#' ## Comparison with v1
#' 
#' Version 1 used \code{signal::resample()} which implements sinc-based

#' bandlimited interpolation. Version 2 uses \code{gsignal::resample()}
#' which implements polyphase filtering (same as MATLAB).
#' 
#' 
#' @examples
#' # Basic usage
#' sig <- sin(seq(0, 4*pi, length.out = 1000))
#' resampled <- resample_signal(sig, target_length = 1200)
#' length(resampled)  # 1200
#' 
#' # With debug output
#' resampled <- resample_signal(sig, target_length = 1050, debug = TRUE)
#' 
#' # Check what method was used
#' attr(resampled, "method_used")
#' attr(resampled, "backend")
#' 
#' # ACIER standard: resample to 18750 samples
#' gait_signal <- rnorm(16000)  # Simulated gait data
#' resampled <- resample_signal(gait_signal, target_length = 18750)
#' 
#' @seealso 
#' \code{\link[gsignal]{resample}}, \code{\link{check_resampling_ratio}}
#' 
#' @export
resample_signal <- function(signal, 
                            target_length,
                            method = c("auto", "direct", "two_stage"),
                            complexity_threshold = 100,
                            intermediate_factor = 1.5,
                            debug = FALSE,
                            use_fallback = TRUE) {
  
  # Match method argument
  method <- match.arg(method)
  
  # -------------------------------------------------------------------------
  # Determine which backend to use
  # -------------------------------------------------------------------------
  backend <- NULL
  
  if (HAS_GSIGNAL_PKG) {
    backend <- "gsignal"
  } else if (HAS_SIGNAL_PKG && use_fallback) {
    backend <- "signal"
    if (debug) {
      message("Note: 'gsignal' not available. Using 'signal' package (less MATLAB-compatible).")
      message("      For best results: install.packages('gsignal')")
    }
  } else if (use_fallback) {
    backend <- "spline"
    if (debug) {
      message("Warning: No signal processing package available.")
      message("         Using spline interpolation (NOT recommended for dynamical analysis).")
      message("         Install gsignal: install.packages('gsignal')")
    }
  } else {
    stop("Package 'gsignal' is required for MATLAB-compatible resampling.\n",
         "Install with: install.packages('gsignal')\n",
         "Or set use_fallback = TRUE to use alternative methods.")
  }
  
  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if (!is.numeric(signal) || length(signal) < 2) {
    stop("'signal' must be a numeric vector with at least 2 elements")
  }
  if (any(!is.finite(signal))) {
    warning("Signal contains non-finite values (NA/NaN/Inf). These may cause issues.")
  }
  if (!is.numeric(target_length) || length(target_length) != 1 || target_length < 1) {
    stop("'target_length' must be a single positive integer")
  }
  if (!is.numeric(complexity_threshold) || complexity_threshold <= 0) {
    stop("'complexity_threshold' must be a positive number")
  }
  if (!is.numeric(intermediate_factor) || intermediate_factor <= 1) {
    stop("'intermediate_factor' must be greater than 1")
  }
  
  original_length <- length(signal)
  target_length <- as.integer(target_length)
  
  # -------------------------------------------------------------------------
  # Edge case: no resampling needed
  # -------------------------------------------------------------------------
  if (original_length == target_length) {
    if (debug) message("No resampling needed: original length equals target length")
    result <- signal
    attr(result, "method_used") <- "none"
    attr(result, "original_length") <- original_length
    attr(result, "filter_phases") <- 1
    attr(result, "was_problematic") <- FALSE
    attr(result, "backend") <- "none"
    return(result)
  }
  
  # -------------------------------------------------------------------------
  # Check if ratio is problematic
  # -------------------------------------------------------------------------
  ratio_check <- check_resampling_ratio(original_length, target_length, complexity_threshold)
  
  if (debug) {
    message(sprintf("Resampling: %d -> %d (ratio = %.4f)",
                    original_length, target_length, ratio_check$ratio))
    message(sprintf("  Backend: %s (polyphase = %s)", 
                    backend, ifelse(backend == "gsignal", "yes", "no")))
    message(sprintf("  GCD = %d, Reduced = %d/%d",
                    ratio_check$gcd, ratio_check$p_reduced, ratio_check$q_reduced))
    message(sprintf("  Filter phases = %d (threshold = %d)",
                    ratio_check$filter_phases, complexity_threshold))
    message(sprintf("  Problematic = %s", ratio_check$is_problematic))
  }
  
  # -------------------------------------------------------------------------
  # Decide on method (direct vs two-stage)
  # -------------------------------------------------------------------------
  use_two_stage <- switch(method,
                          "auto" = ratio_check$is_problematic,
                          "two_stage" = TRUE,
                          "direct" = FALSE)
  
  # -------------------------------------------------------------------------
  # Helper: check if gsignal::resample(x, p, q) would overflow R's 32-bit
  # integer limit, causing "cannot allocate vector of negative length"
  #
  # ROOT CAUSE :
  #   gsignal::resample does NOT reduce p/q by their GCD before calling
  #   upfirdn(). Internally, upfirdn allocates a buffer of size:
  #       (length(x) - 1) * p + length(h)
  #   where h is the FIR filter (length ~ 2*10*max(p,q) + 1).
  #   R vectors are limited to 2^31 - 1 elements (32-bit indexing).
  #   When (N-1)*p exceeds this, the integer wraps to negative, and R
  #   attempts to allocate a negative-length vector -> crash.
  #
  # EXAMPLE :
  #   SF = 1.668 Hz -> truncation_index = 38370
  #   Two-stage Stage 1: resample(x, p=57555, q=38370)
  #     upfirdn buffer = (38370-1)*57555 = 2,209 M > 2^31 = 2,147 M -> CRASH
  #   Direct: resample(x, p=18750, q=38370)
  #     upfirdn buffer = (38370-1)*18750 = 719 M < 2^31 -> OK
  #
  # The two-stage workaround creates a LARGER p than direct, making
  # overflow MORE likely. The fix detects this and falls back to direct.
  # -------------------------------------------------------------------------
  .would_overflow <- function(sig_len, p, q) {
    # Check unreduced p (as gsignal actually passes to upfirdn)
    # Conservative: check (sig_len - 1) * p (dominant term; filter length
    # adds ~max(p,q)*20 which is negligible relative to the product)
    return(as.double(sig_len - 1) * as.double(p) > .Machine$integer.max)
  }
  
  # -------------------------------------------------------------------------
  # Helper function for actual resampling
  # -------------------------------------------------------------------------
  do_resample <- function(x, target_len, orig_len = length(x)) {
    result <- switch(backend,
                     "gsignal" = .resample_polyphase(x, target_len, orig_len),
                     "signal" = .resample_bandlimited(x, target_len, orig_len),
                     "spline" = .resample_spline(x, target_len))
    return(as.numeric(result))
  }
  
  # -------------------------------------------------------------------------
  # Safety: if two-stage would cause integer overflow, try alternatives
  # This is the key fix for cases where
  # the "problematic ratio" workaround creates worse overflow than direct.
  # -------------------------------------------------------------------------
  if (use_two_stage && backend %in% c("gsignal", "signal")) {
    intermediate_length_check <- as.integer(ceiling(original_length * intermediate_factor))
    if (.would_overflow(original_length, intermediate_length_check, original_length)) {
      # Two-stage Stage 1 would overflow -> try direct instead
      if (!.would_overflow(original_length, target_length, original_length)) {
        if (debug) {
          ts_buf <- as.double(original_length - 1) * as.double(intermediate_length_check)
          dr_buf <- as.double(original_length - 1) * as.double(target_length)
          message(sprintf("  [!] Two-stage Stage 1 would overflow R integer limit:"))
          message(sprintf("      upfirdn buffer = (N-1)*p = %.0f > 2^31 = %.0f", 
                          ts_buf, .Machine$integer.max))
          message(sprintf("      Falling back to direct resampling (buffer = %.0f, safe)", dr_buf))
        }
        use_two_stage <- FALSE
      } else {
        # Both would overflow -> fall back to spline
        if (debug) {
          message("  [!] Both two-stage and direct would overflow R integer limit")
          message("      Falling back to spline interpolation")
        }
        backend <- "spline"
        use_two_stage <- FALSE
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # Perform resampling (with tryCatch for safety)
  # -------------------------------------------------------------------------
  if (use_two_stage) {
    if (debug) message("  Using two-stage resampling")
    
    # Stage 1: Upsample to intermediate length
    intermediate_length <- as.integer(ceiling(original_length * intermediate_factor))
    
    if (debug) {
      message(sprintf("  Stage 1: %d -> %d (factor = %.2f)",
                      original_length, intermediate_length, intermediate_factor))
    }
    
    result <- tryCatch({
      temp <- do_resample(signal, intermediate_length, original_length)
      
      # Stage 2: Resample to target
      if (debug) {
        diag_stage2 <- check_resampling_ratio(length(temp), target_length, complexity_threshold)
        message(sprintf("  Stage 2: %d -> %d (phases = %d)",
                        length(temp), target_length, diag_stage2$filter_phases))
      }
      
      do_resample(temp, target_length, length(temp))
    }, error = function(e) {
      # Two-stage failed (e.g. memory/overflow) -> try direct, then spline
      if (debug) message(sprintf("  Two-stage failed (%s); trying direct...", e$message))
      tryCatch(
        do_resample(signal, target_length, original_length),
        error = function(e2) {
          if (debug) message(sprintf("  Direct also failed (%s); using spline fallback", e2$message))
          .resample_spline(signal, target_length)
        }
      )
    })
    method_used <- "two_stage"
    
  } else {
    if (debug) message("  Using direct resampling")
    result <- tryCatch(
      do_resample(signal, target_length, original_length),
      error = function(e) {
        if (debug) message(sprintf("  Direct failed (%s); using spline fallback", e$message))
        .resample_spline(signal, target_length)
      }
    )
    method_used <- "direct"
  }
  
  # -------------------------------------------------------------------------
  # Ensure exact length (handle potential rounding errors)
  # -------------------------------------------------------------------------
  if (length(result) != target_length) {
    if (debug) {
      message(sprintf("  Length adjustment: %d -> %d", length(result), target_length))
    }
    if (length(result) > target_length) {
      result <- result[1:target_length]
    } else {
      # Pad with last value if too short (rare edge case)
      result <- c(result, rep(tail(result, 1), target_length - length(result)))
    }
  }
  
  # -------------------------------------------------------------------------
  # Attach metadata as attributes
  # -------------------------------------------------------------------------
  attr(result, "method_used") <- method_used
  attr(result, "original_length") <- original_length
  attr(result, "filter_phases") <- ratio_check$filter_phases
  attr(result, "was_problematic") <- ratio_check$is_problematic
  attr(result, "backend") <- backend
  
  return(result)
}

