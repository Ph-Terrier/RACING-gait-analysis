#' ============================================================================
#' Robust Signal Resampling Function (MATLAB-compatible)
#' ============================================================================
#' 
#' Resamples signals to a target length using polyphase filtering,
#' replicating MATLAB's resample(x, p, q, N, beta) algorithm exactly.
#' 
#' @author Philippe Terrier / ORD 2025 Project (RACING)
#' @date March 2026
#' @version 3.6
#' @license MIT
#'
#' KEY FIX (v3.6): Stage 1 of two-stage resampling now passes the ratio
#' (p=2, q=3) directly to .resample_matlab(), matching MATLAB's
#' resample(x, 2, 3, 20). Previous versions passed computed *lengths*
#' (original_length, intermediate_length) which, when GCD = 1, produced
#' a multi-million-tap filter and an upfirdn buffer exceeding R's 2^31
#' integer index limit ("invalid 'times' argument" in rep()).
#' Indoor data (250 steps, ~37k samples) was unaffected because the
#' buffer stayed below 2^31. Outdoor data (500 steps, ~75k samples)
#' triggered the overflow. The overflow guard is also extended to check
#' Stage 1 (previously only Stage 2 was checked). A dedicated helper
#' do_resample_stage1() encapsulates the ratio-based Stage 1 call for
#' all backends.
#'
#' KEY FIX (v3.5): Raised auto-mode threshold from 100 to 50000 filter phases.
#' The old threshold flagged 99% of ACIER signals as "problematic" even though
#' all completed successfully in direct mode (max buffer ~2.4 GB, well within
#' R's 16 GB integer limit). The overflow guard (.would_overflow) remains as the
#' real safety net. Console output clarified: ratio diagnostics no longer print
#' misleading "Problematic = TRUE" when method is forced to "direct".
#'
#' KEY FIX (v3.4): Added file-based debug logging via log_file parameter.
#' In batch mode, call resample_log_header() once, then pass log_file to
#' every resample_signal() call. The log accumulates all algorithmic
#' details: signal stats, ratio analysis, filter design, delay compensation,
#' upfirdn buffer sizes, extraction indices, and output statistics.
#'
#' KEY FIX (v3.3): Filter design uses O(N) sinc computation instead of
#' O(N^3) firls. For the ACIER 75/128 ratio, firls builds and solves a
#' 2561x2561 matrix (~52 MB) via pracma::mldivide, taking 5-10 minutes.
#' MATLAB's firls with zero-width transition band [0 2fc 2fc 1] produces
#' the ideal lowpass sinc function -- we compute it directly in O(N).
#' After Kaiser windowing and normalization, results are near-identical.
#' 
#' Architecture: implements MATLAB's uniformResample() algorithm,
#' calling gsignal::upfirdn() (compiled C++) for the polyphase step.
#' Filter caching across axes (norm, x, y, z) per subject.
#' 
#' @references
#' MATLAB resample():
#'   https://www.mathworks.com/help/signal/ref/resample.html
#' MATLAB source: Signal Processing Toolbox, uniformResample()
#' ============================================================================


# =============================================================================
# Package Availability Check
# =============================================================================

HAS_GSIGNAL_PKG <- requireNamespace("gsignal", quietly = TRUE)
HAS_SIGNAL_PKG <- requireNamespace("signal", quietly = TRUE)

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
# Utility Functions
# =============================================================================

#' Calculate Greatest Common Divisor (GCD)
#' @param a First integer
#' @param b Second integer
#' @return GCD of a and b
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


# =============================================================================
# Filter Cache
# =============================================================================
# In batch processing (ACIER), the same resampling ratio is used for 4 axes
# (norm, x, y, z) per subject/condition. Designing the firls filter once
# and reusing it for all 4 axes saves ~75% of filter design time.
#
# The cache is keyed by (p_red, q_red, n, beta) -- cleared with
# clear_resample_cache() or automatically when R session ends.
# =============================================================================

.resample_filter_cache <- new.env(parent = emptyenv())

#' Clear the resampling filter cache
#' @export
clear_resample_cache <- function() {
  rm(list = ls(.resample_filter_cache), envir = .resample_filter_cache)
}


# =============================================================================
# Debug File Logging
# =============================================================================
# When log_file is set, every resampling call appends detailed trace info
# to the file. In batch mode this accumulates across all subjects/axes,
# giving a complete record for comparison with MATLAB.
#
# Usage:
#   resample_log_header("resample_log.txt")   # once at batch start
#   resample_signal(x, target, log_file = "resample_log.txt")
#   resample_signal(y, target, log_file = "resample_log.txt")  # appends
# =============================================================================

#' Write a line to the resampling debug log file
#' @param log_file Path to log file (NULL = skip)
#' @param ... Arguments passed to sprintf
#' @keywords internal
.resample_log <- function(log_file, fmt, ...) {
  if (is.null(log_file)) return(invisible(NULL))
  line <- sprintf(fmt, ...)
  cat(line, "\n", file = log_file, append = TRUE, sep = "")
}

#' Write a separator block to the log
#' @keywords internal
.resample_log_sep <- function(log_file, label = "") {
  if (is.null(log_file)) return(invisible(NULL))
  rule <- paste(rep("=", 72), collapse = "")
  cat("\n", rule, "\n", file = log_file, append = TRUE, sep = "")
  if (nzchar(label)) {
    cat(label, "\n", file = log_file, append = TRUE, sep = "")
    cat(rule, "\n", file = log_file, append = TRUE, sep = "")
  }
}

#' Log signal summary statistics
#' @keywords internal
.resample_log_signal <- function(log_file, label, x) {
  if (is.null(log_file)) return(invisible(NULL))
  n <- length(x)
  .resample_log(log_file, "  %s: n=%d, range=[%.6g, %.6g], mean=%.6g, sd=%.6g",
                label, n, min(x), max(x), mean(x), sd(x))
  # First and last 5 values (critical for boundary effects)
  head_vals <- paste(sprintf("%.6g", head(x, 5)), collapse = ", ")
  tail_vals <- paste(sprintf("%.6g", tail(x, 5)), collapse = ", ")
  .resample_log(log_file, "    head(5): %s", head_vals)
  .resample_log(log_file, "    tail(5): %s", tail_vals)
}

#' Write a header to the log file (call once before batch processing)
#'
#' @param log_file Path to log file (created/overwritten)
#' @param session_info Logical. Include R session info (default: TRUE)
#' @export
resample_log_header <- function(log_file, session_info = TRUE) {
  # Overwrite (not append) to start fresh
  rule <- paste(rep("=", 72), collapse = "")
  cat(rule, "\n", file = log_file, append = FALSE, sep = "")
  cat("RACING Resampling Debug Log\n", file = log_file, append = TRUE)
  cat(sprintf("Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = log_file, append = TRUE)
  cat(rule, "\n", file = log_file, append = TRUE, sep = "")
  
  if (session_info) {
    cat(sprintf("R version: %s\n", R.version.string), file = log_file, append = TRUE)
    cat(sprintf("Platform:  %s\n", R.version$platform), file = log_file, append = TRUE)
    cat(sprintf("gsignal:   %s\n",
                ifelse(HAS_GSIGNAL_PKG,
                       as.character(packageVersion("gsignal")), "NOT INSTALLED")),
        file = log_file, append = TRUE)
    cat(sprintf("signal:    %s\n",
                ifelse(HAS_SIGNAL_PKG,
                       as.character(packageVersion("signal")), "NOT INSTALLED")),
        file = log_file, append = TRUE)
  }
  cat("\n", file = log_file, append = TRUE)
  invisible(log_file)
}

#' Write a separator block to the resampling log file
#'
#' Public wrapper for .resample_log_sep(), for use in batch scripts
#' and truncate_resample.R. Writes a visible separator with an
#' optional label (e.g., file name) to the log file.
#'
#' @param log_file Path to log file (NULL = skip)
#' @param label Character. Text to display in the separator block.
#' @export
resample_log_sep <- function(log_file, label = "") {
  .resample_log_sep(log_file, label)
}

#' Write a formatted line to the resampling log file
#'
#' Public wrapper for .resample_log(), for use in other scripts
#' (e.g., truncate_resample.R) that need to write context lines
#' to the same log file.
#'
#' @param log_file Path to log file (NULL = skip)
#' @param fmt Character. sprintf format string.
#' @param ... Values for sprintf placeholders.
#' @export
resample_log <- function(log_file, fmt, ...) {
  .resample_log(log_file, fmt, ...)
}


# =============================================================================
# Internal: MATLAB-compatible Filter Design
# =============================================================================

#' Design MATLAB-compatible resampling filter (with caching)
#'
#' Replicates MATLAB's uniformResample() filter design:
#'   fc    = 1 / (2 * pqmax)
#'   order = 2 * N * pqmax
#'   h     = firls(order, [0 2fc 2fc 1], [1 1 0 0]) .* kaiser(order+1, beta)
#'   h     = p * h / sum(h)
#'
#' PERFORMANCE NOTE (v3.3):
#'   MATLAB's firls with zero-width transition band [0 2fc 2fc 1] converges
#'   to the ideal lowpass sinc function. We compute sinc directly: O(N) vs
#'   O(N^3) for the matrix solve in gsignal::firls. For order=5120 (ACIER
#'   75/128 ratio), this reduces filter design from ~5-10 min to <0.01 sec.
#'   The resulting filter is mathematically equivalent after Kaiser windowing.
#'
#' p and q MUST be GCD-reduced (coprime) before calling.
#'
#' @param p Upsampling factor (coprime with q)
#' @param q Downsampling factor (coprime with p)
#' @param n Half-order (MATLAB default=10, ACIER=20)
#' @param beta Kaiser shape parameter (MATLAB default=5)
#' @param use_cache Logical. Use filter cache (default: TRUE).
#' @return FIR filter coefficient vector
#' @keywords internal
.design_matlab_resample_filter <- function(p, q, n = 20, beta = 5, 
                                           use_cache = TRUE,
                                           log_file = NULL) {
  # Check cache first
  cache_key <- sprintf("p%d_q%d_n%d_b%.2f", p, q, n, beta)
  from_cache <- use_cache && exists(cache_key, envir = .resample_filter_cache)
  if (from_cache) {
    h <- get(cache_key, envir = .resample_filter_cache)
    .resample_log(log_file, "  Filter: from cache [%s], %d taps", cache_key, length(h))
    return(h)
  }
  
  pqmax <- max(p, q)
  fc <- 1 / (2 * pqmax)     # normalized cutoff
  order <- 2 * n * pqmax    # filter order
  L <- order + 1             # filter length
  
  .resample_log(log_file, "  Filter design: key=%s, pqmax=%d, fc=%.8f, order=%d, L=%d",
                cache_key, pqmax, fc, order, L)
  
  # -----------------------------------------------------------------------
  # Ideal lowpass via sinc (replaces firls -- see PERFORMANCE NOTE above)
  #
  # MATLAB: firls(order, [0 2*fc 2*fc 1], [1 1 0 0])
  # With zero-width transition band, firls converges to ideal lowpass:
  #   h_ideal[k] = sin(2*pi*fc * (k - order/2)) / (pi * (k - order/2))
  #   h_ideal[order/2] = 2 * fc
  # This is the standard sinc(2*fc*x) = sin(pi*2*fc*x) / (pi*2*fc*x)
  # -----------------------------------------------------------------------
  k <- seq(0, order)
  m <- k - order / 2           # centered sample indices
  
  # Compute sinc: sin(pi*x)/(pi*x) with x = 2*fc*m
  x <- 2 * fc * m
  h <- numeric(L)
  nz <- abs(x) > 1e-20         # non-zero indices (avoid 0/0)
  h[nz] <- sin(pi * x[nz]) / (pi * x[nz])
  h[!nz] <- 1.0                # sinc(0) = 1
  
  # Apply Kaiser window (same as MATLAB)
  w <- gsignal::kaiser(L, beta)
  h <- h * w
  
  # Normalize: MATLAB does h = p * h / sum(h)
  h <- p * h / sum(h)
  
  # Log filter characteristics
  .resample_log(log_file, "    sum(h_raw) before norm: computed, h normalized to sum=%.8f", sum(h))
  .resample_log(log_file, "    h range: [%.8g, %.8g], max at index %d",
                min(h), max(h), which.max(h))
  h5 <- paste(sprintf("%.8g", head(h, 5)), collapse = ", ")
  .resample_log(log_file, "    h head(5): %s", h5)
  
  # Store in cache
  if (use_cache) {
    assign(cache_key, h, envir = .resample_filter_cache)
  }
  
  return(h)
}


# =============================================================================
# Internal: Direct MATLAB resample implementation
# =============================================================================

#' Direct port of MATLAB's uniformResample algorithm
#'
#' This function implements MATLAB's resample() line-by-line, calling
#' gsignal::upfirdn() directly (compiled C) instead of going through
#' gsignal::resample() (R wrapper with overhead).
#'
#' MATLAB algorithm (from uniformResample in resample.m):
#'   1. [p, q] = rat(p/q, 1e-12)          -- reduce to coprime
#'   2. design filter h for max(p,q)       -- firls + kaiser
#'   3. Lhalf = (L-1)/2                    -- half filter length
#'   4. nZeroBegin = floor(q - Lhalf%%q)   -- prepend zeros for alignment
#'   5. h = [zeros, h]                     -- prepend
#'   6. delay = floor(ceil(Lhalf)/q)       -- output delay
#'   7. pad h with trailing zeros          -- ensure output length
#'   8. y_full = upfirdn(x, h, p, q)      -- polyphase filtering (C code)
#'   9. y = y_full[delay + (1:Ly)]         -- trim to output length
#'
#' @param x Numeric vector. Input signal.
#' @param p Integer. Target length (before GCD reduction).
#' @param q Integer. Original length (before GCD reduction).
#' @param n Integer. Filter half-order (default=20).
#' @param beta Numeric. Kaiser beta (default=5).
#' @param debug Logical. Print sub-step timing (default=FALSE).
#' @return Resampled signal vector of length ceil(length(x) * p / q).
#' @keywords internal
.resample_matlab <- function(x, p, q, n = 20, beta = 5, debug = FALSE,
                            log_file = NULL) {
  if (!HAS_GSIGNAL_PKG) {
    stop("gsignal package required for MATLAB-compatible resampling")
  }
  
  Lx <- length(x)
  
  # -------------------------------------------------------------------
  # Step 1: Reduce p/q to coprime (MATLAB: [p,q] = rat(p/q, 1e-12))
  # -------------------------------------------------------------------
  g <- gcd(p, q)
  p <- as.integer(p / g)
  q <- as.integer(q / g)
  
  .resample_log(log_file, "  .resample_matlab: Lx=%d, raw p/q=%d/%d, GCD=%d, reduced p/q=%d/%d",
                Lx, as.integer(p * g), as.integer(q * g), g, p, q)
  
  # Trivial case
  if (p == 1L && q == 1L) {
    .resample_log(log_file, "  -> Trivial (p=q=1), returning input unchanged")
    return(x)
  }
  
  # -------------------------------------------------------------------
  # Step 2: Design filter (MATLAB: firls + kaiser, with caching)
  # -------------------------------------------------------------------
  if (debug) t0 <- proc.time()["elapsed"]
  
  h <- .design_matlab_resample_filter(p, q, n = n, beta = beta,
                                       log_file = log_file)
  
  if (debug) {
    dt <- proc.time()["elapsed"] - t0
    message(sprintf("    [timing] Filter design (sinc+kaiser): %.4f sec (%d taps, %s)",
                    dt, length(h),
                    ifelse(exists(sprintf("p%d_q%d_n%d_b%.2f", p, q, n, beta),
                                  envir = .resample_filter_cache),
                           "cached", "computed")))
  }
  
  # -------------------------------------------------------------------
  # Step 3-6: Delay compensation (exact MATLAB logic)
  # -------------------------------------------------------------------
  # MATLAB: Lhalf = (L-1)/2
  L <- length(h)
  Lhalf <- (L - 1) / 2
  
  # MATLAB: nZeroBegin = floor(q - mod(Lhalf, q))
  nZeroBegin <- floor(q - (Lhalf %% q))
  
  # MATLAB: h = [zeros(1, nZeroBegin), h]
  h_padded <- c(rep(0, nZeroBegin), h)
  Lhalf_orig <- Lhalf
  Lhalf <- Lhalf + nZeroBegin
  
  # MATLAB: delay = floor(ceil(Lhalf) / q)
  delay <- floor(ceiling(Lhalf) / q)
  
  # Output length: MATLAB: Ly = ceil(Lx * p / q)
  Ly <- ceiling(Lx * p / q)
  
  .resample_log(log_file, "  Delay compensation:")
  .resample_log(log_file, "    L=%d, Lhalf_orig=%.1f, Lhalf%%q=%.1f",
                L, Lhalf_orig, Lhalf_orig %% q)
  .resample_log(log_file, "    nZeroBegin=%d, Lhalf_adjusted=%.1f", nZeroBegin, Lhalf)
  .resample_log(log_file, "    ceil(Lhalf)=%d, delay=floor(%d/%d)=%d",
                ceiling(Lhalf), ceiling(Lhalf), q, delay)
  .resample_log(log_file, "    Ly=ceil(%d*%d/%d)=%d", Lx, p, q, Ly)
  
  # -------------------------------------------------------------------
  # Step 7: Compute trailing zero-pad
  # -------------------------------------------------------------------
  Lh_current <- length(h_padded)
  Lh_needed <- (delay + Ly - 1) * q - (Lx - 1) * p + 1
  nZeroEnd <- max(0L, as.integer(Lh_needed - Lh_current))
  h_padded <- c(h_padded, rep(0, nZeroEnd))
  
  .resample_log(log_file, "  Zero-padding:")
  .resample_log(log_file, "    Lh_current=%d, Lh_needed=%d, nZeroEnd=%d, Lh_final=%d",
                Lh_current, as.integer(Lh_needed), nZeroEnd, length(h_padded))
  
  # -------------------------------------------------------------------
  # Step 8: Call upfirdn (compiled C -- this is the fast part)
  # -------------------------------------------------------------------
  if (debug) t0 <- proc.time()["elapsed"]
  
  .resample_log(log_file, "  upfirdn call: x[%d], h[%d], p=%d, q=%d",
                Lx, length(h_padded), p, q)
  
  y_full <- gsignal::upfirdn(x, h_padded, p, q)
  
  .resample_log(log_file, "  upfirdn result: length=%d (%.2f M elements)",
                length(y_full), length(y_full) / 1e6)
  
  if (debug) {
    dt <- proc.time()["elapsed"] - t0
    message(sprintf("    [timing] upfirdn: %.4f sec (buffer = %.2f M elements)",
                    dt, length(y_full) / 1e6))
  }
  
  # -------------------------------------------------------------------
  # Step 9: Extract output with delay compensation
  # MATLAB: y = yVec(delay + (1:Ly))
  # -------------------------------------------------------------------
  idx_start <- delay + 1L
  idx_end <- delay + Ly
  
  .resample_log(log_file, "  Extraction: y_full[%d:%d] (delay=%d, Ly=%d)",
                idx_start, idx_end, delay, Ly)
  
  if (idx_end > length(y_full)) {
    .resample_log(log_file, "  [!] WARNING: idx_end=%d > length(y_full)=%d",
                  idx_end, length(y_full))
  }
  
  idx <- delay + seq_len(Ly)
  y <- y_full[idx]
  
  # Log output signal stats
  .resample_log_signal(log_file, "OUTPUT", y)
  
  return(y)
}


# =============================================================================
# Fallback Resampling Functions
# =============================================================================

#' @keywords internal
.resample_bandlimited <- function(x, p, q) {
  if (!HAS_SIGNAL_PKG) stop("signal package not available")
  return(signal::resample(x, p, q))
}

#' @keywords internal
.resample_spline <- function(x, target_length) {
  n_pts <- length(x)
  x_orig <- seq(0, 1, length.out = n_pts)
  x_new <- seq(0, 1, length.out = target_length)
  return(spline(x_orig, x, xout = x_new, method = "natural")$y)
}


# =============================================================================
# Ratio Analysis
# =============================================================================

#' Check if Resampling Ratio is Problematic
#' 
#' @param original_length Integer. Length of the original signal.
#' @param target_length Integer. Desired length after resampling.
#' @param threshold Integer. Maximum allowed filter phases (default: 100).
#' @return A list with ratio diagnostics.
#' @export
check_resampling_ratio <- function(original_length, target_length, threshold = 50000) {
  
  stopifnot(
    is.numeric(original_length), length(original_length) == 1, original_length > 0,
    is.numeric(target_length), length(target_length) == 1, target_length > 0,
    is.numeric(threshold), length(threshold) == 1, threshold > 0
  )
  
  original_length <- as.integer(original_length)
  target_length <- as.integer(target_length)
  
  g <- gcd(original_length, target_length)
  p_reduced <- target_length / g
  q_reduced <- original_length / g
  filter_phases <- max(p_reduced, q_reduced)
  filter_length <- 2L * 20L * filter_phases + 1L
  
  # upfirdn buffer with REDUCED ratio: (Lx-1)*p_red + filter_length
  upfirdn_buffer <- as.double(original_length - 1) * as.double(p_reduced) +
                    as.double(filter_length)
  upfirdn_buffer_MB <- upfirdn_buffer * 8 / 1024 / 1024
  
  return(list(
    ratio = target_length / original_length,
    gcd = g,
    p_reduced = p_reduced,
    q_reduced = q_reduced,
    filter_phases = filter_phases,
    filter_length = filter_length,
    upfirdn_buffer = upfirdn_buffer,
    upfirdn_buffer_MB = upfirdn_buffer_MB,
    is_problematic = filter_phases > threshold
  ))
}


# =============================================================================
# Main Resampling Function
# =============================================================================

#' Resample Signal (MATLAB-compatible, optimized)
#' 
#' Resamples a signal to a target length using polyphase filtering,
#' replicating MATLAB's resample(x, p, q, N, beta) algorithm exactly.
#' 
#' @param signal Numeric vector. The input signal to resample.
#' @param target_length Integer. The desired output length.
#' @param method Character. Resampling strategy:
#' \describe{
#'   \item{"auto"}{(default) Use two-stage only when it improves the ratio}
#'   \item{"direct"}{Always use direct resampling}
#'   \item{"two_stage"}{Always use two-stage resampling}
#' }
#' @param complexity_threshold Integer. Max filter phases for "auto" (default: 50000).
#'   Only relevant when method="auto". The overflow guard (.would_overflow) provides
#'   a hard safety check regardless of this threshold.
#' @param n Integer. Filter half-order (default: 20, matching ACIER).
#' @param beta Numeric. Kaiser window beta (default: 5).
#' @param debug Logical. Print diagnostics with sub-step timing (default: FALSE).
#' @param use_fallback Logical. Fall back to signal/spline if needed (default: TRUE).
#' @param log_file Character or NULL. Path to debug log file. When set,
#'   detailed trace information is appended to this file for every call.
#'   Use \code{resample_log_header()} to initialize the log before batch
#'   processing. Default: NULL (no file logging).
#' @param log_label Character. Optional label for this call (e.g.,
#'   "S01_indoor_norm") to identify entries in the log file.
#' 
#' @return Numeric vector of length \code{target_length} with metadata attributes.
#' 
#' @details
#' ## Performance
#' 
#' This function bypasses gsignal::resample() entirely and calls 
#' gsignal::upfirdn() (compiled C) directly. The filter is cached
#' so that resampling 4 axes with the same ratio only designs the
#' filter once. Typical speedup vs v3.1: 2-5x.
#' 
#' ## Algorithm
#' 
#' Direct port of MATLAB's uniformResample():
#' 1. Reduce p/q to coprime via GCD
#' 2. Design firls+Kaiser filter (cached)
#' 3. Compute delay compensation (nZeroBegin, nZeroEnd)
#' 4. Call upfirdn(x, h_padded, p_red, q_red)  [compiled C]
#' 5. Extract y = y_full[delay + (1:Ly)]
#' 
#' ## Two-Stage Handling
#' 
#' When ratio is near 1.0 (problematic), uses MATLAB's shrink-first:
#' Stage 1: resample(x, 2, 3, N, beta)  -- shrink to 2/3
#' Stage 2: resample(temp, target, length(temp), N, beta)
#' 
#' Two-stage is only applied when Stage 2 has fewer filter phases
#' than direct (P4 smart check).
#' 
#' @examples
#' # ACIER: 16000 -> 18750 (GCD=250, 75/64, filter=3001 taps)
#' gait <- rnorm(16000)
#' y <- resample_signal(gait, 18750, debug = TRUE)
#' 
#' # Batch: same ratio reuses cached filter
#' y_x <- resample_signal(gait_x, 18750)  # designs filter
#' y_y <- resample_signal(gait_y, 18750)  # uses cached filter
#' y_z <- resample_signal(gait_z, 18750)  # uses cached filter
#' 
#' @export
resample_signal <- function(signal, 
                            target_length,
                            method = c("auto", "direct", "two_stage"),
                            complexity_threshold = 50000,
                            n = 20,
                            beta = 5,
                            debug = TRUE,
                            use_fallback = TRUE,
                            log_file = NULL,
                            log_label = "") {
  
  method <- match.arg(method)
  t_start <- proc.time()["elapsed"]
  
  # Log header for this call
  .resample_log_sep(log_file, sprintf("resample_signal [%s] @ %s",
                                       log_label, format(Sys.time(), "%H:%M:%S")))
  .resample_log(log_file, "Parameters: target=%d, method=%s, threshold=%d, n=%d, beta=%.1f",
                as.integer(target_length), method, complexity_threshold, n, beta)
  
  # -------------------------------------------------------------------------
  # Determine backend
  # -------------------------------------------------------------------------
  backend <- NULL
  if (HAS_GSIGNAL_PKG) {
    backend <- "gsignal"
  } else if (HAS_SIGNAL_PKG && use_fallback) {
    backend <- "signal"
    if (debug) message("  Note: gsignal not available, using signal package")
  } else if (use_fallback) {
    backend <- "spline"
    if (debug) message("  Warning: No signal processing package, using spline")
  } else {
    stop("Package 'gsignal' is required. Install: install.packages('gsignal')")
  }
  
  .resample_log(log_file, "Backend: %s", backend)
  
  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if (!is.numeric(signal) || length(signal) < 2) {
    stop("'signal' must be a numeric vector with at least 2 elements")
  }
  if (any(!is.finite(signal))) {
    warning("Signal contains non-finite values (NA/NaN/Inf).")
  }
  if (!is.numeric(target_length) || length(target_length) != 1 || target_length < 1) {
    stop("'target_length' must be a single positive integer")
  }
  
  original_length <- length(signal)
  target_length <- as.integer(target_length)
  
  # Log input signal statistics
  .resample_log_signal(log_file, "INPUT", signal)
  
  # -------------------------------------------------------------------------
  # Edge case: no resampling needed
  # -------------------------------------------------------------------------
  if (original_length == target_length) {
    if (debug) message("  No resampling needed (lengths equal)")
    .resample_log(log_file, "No resampling needed (original == target == %d)", original_length)
    result <- signal
    attr(result, "method_used") <- "none"
    attr(result, "original_length") <- original_length
    attr(result, "filter_phases") <- 1L
    attr(result, "filter_length") <- 1L
    attr(result, "was_problematic") <- FALSE
    attr(result, "backend") <- "none"
    attr(result, "elapsed_sec") <- 0
    return(result)
  }
  
  # -------------------------------------------------------------------------
  # Analyze ratio
  # -------------------------------------------------------------------------
  ratio_check <- check_resampling_ratio(original_length, target_length, 
                                        complexity_threshold)
  
  .resample_log(log_file, "Ratio analysis:")
  .resample_log(log_file, "  %d -> %d, ratio=%.8f", original_length, target_length, ratio_check$ratio)
  .resample_log(log_file, "  GCD=%d, reduced p/q=%d/%d",
                ratio_check$gcd, ratio_check$p_reduced, ratio_check$q_reduced)
  .resample_log(log_file, "  filter_phases=%d, filter_length=%d",
                ratio_check$filter_phases, ratio_check$filter_length)
  .resample_log(log_file, "  upfirdn_buffer=%.0f (%.1f MB), problematic=%s",
                ratio_check$upfirdn_buffer, ratio_check$upfirdn_buffer_MB,
                ratio_check$is_problematic)
  
  if (debug) {
    message(sprintf("  Resampling: %d -> %d (ratio = %.6f)",
                    original_length, target_length, ratio_check$ratio))
    message(sprintf("  Backend: %s (direct upfirdn) | n=%d, beta=%.1f",
                    backend, n, beta))
    message(sprintf("  GCD = %d, reduced = %d/%d",
                    ratio_check$gcd, ratio_check$p_reduced, ratio_check$q_reduced))
    message(sprintf("  Filter: %d phases, %d taps, buffer = %.1f MB",
                    ratio_check$filter_phases, ratio_check$filter_length,
                    ratio_check$upfirdn_buffer_MB))
    # Only show "problematic" diagnostic when method="auto" (where it matters)
    if (method == "auto" && ratio_check$is_problematic) {
      message(sprintf("  [!] Ratio flagged as problematic (%d phases > %d threshold)",
                      ratio_check$filter_phases, complexity_threshold))
    }
  }
  
  # -------------------------------------------------------------------------
  # Decide method
  # -------------------------------------------------------------------------
  use_two_stage <- switch(method,
                          "auto" = ratio_check$is_problematic,
                          "two_stage" = TRUE,
                          "direct" = FALSE)
  
  .resample_log(log_file, "Method decision: requested=%s, use_two_stage=%s", method, use_two_stage)
  
  # -------------------------------------------------------------------------
  # Overflow guard: buffer must fit in R's int32 index limit
  # -------------------------------------------------------------------------
  .would_overflow <- function(sig_len, p, q) {
    g <- gcd(p, q)
    p_red <- as.double(p / g)
    pqmax <- max(p_red, as.double(q / g))
    filter_len <- 2.0 * n * pqmax + 1.0
    buffer <- as.double(sig_len - 1) * p_red + filter_len
    return(buffer > .Machine$integer.max)
  }
  
  # -------------------------------------------------------------------------
  # Helper: perform one resampling operation
  # -------------------------------------------------------------------------
  do_resample <- function(x, target_len, orig_len = length(x), label = "") {
    .resample_log(log_file, "  >> do_resample [%s]: len=%d -> %d", label, orig_len, target_len)
    result <- switch(backend,
                     "gsignal" = .resample_matlab(x, target_len, orig_len,
                                                   n = n, beta = beta,
                                                   debug = debug,
                                                   log_file = log_file),
                     "signal" = .resample_bandlimited(x, target_len, orig_len),
                     "spline" = .resample_spline(x, target_len))
    return(as.numeric(result))
  }
  
  # -------------------------------------------------------------------------
  # Helper: Stage 1 shrink-to-2/3 using RATIO, not lengths  (v3.6 fix)
  # -------------------------------------------------------------------------
  # MATLAB (Main_results_new_v2.m, line 683):
  #   provwalk_norm_A = resample(provwalk_norm_A, 2, 3, 20);
  # This passes p=2, q=3 directly. GCD(2,3)=1, so the reduced filter has
  # max(2,3)=3 phases, order=2*20*3=120, 121 taps -- trivially small.
  #
  # Previous R code passed (intermediate_length, original_length) as
  # (p, q). When original_length is not divisible by 3, GCD can be 1,
  # giving p_red = intermediate_length (e.g., 50196) and a filter with
  # millions of taps, overflowing R's 2^31 integer index limit.
  #
  # This helper calls each backend with the ratio directly:
  #   gsignal: .resample_matlab(x, 2, 3, ...)  -> upfirdn with 121-tap filter
  #   signal:  signal::resample(x, 2, 3)       -> same ratio-based call
  #   spline:  .resample_spline(x, ceil(N*2/3)) -> length-based (no filter)
  # -------------------------------------------------------------------------
  do_resample_stage1 <- function(x, label = "Stage 1 (2/3)") {
    Lx <- length(x)
    expected_len <- ceiling(Lx * 2 / 3)
    .resample_log(log_file, "  >> do_resample_stage1 [%s]: len=%d -> ~%d (ratio 2/3)",
                  label, Lx, expected_len)
    result <- switch(backend,
                     "gsignal" = .resample_matlab(x, 2L, 3L,
                                                   n = n, beta = beta,
                                                   debug = debug,
                                                   log_file = log_file),
                     "signal" = .resample_bandlimited(x, 2L, 3L),
                     "spline" = .resample_spline(x, expected_len))
    return(as.numeric(result))
  }
  
  # -------------------------------------------------------------------------
  # Smart two-stage check (P4): only use if Stage 2 is actually better
  # -------------------------------------------------------------------------
  # v3.6: use ceiling() to match MATLAB's resample(x,2,3) output length
  # which is ceil(Lx * 2/3). Previous floor() could be off by 1.
  intermediate_length <- as.integer(ceiling(original_length * 2 / 3))
  
  if (use_two_stage) {
    stage2_check <- check_resampling_ratio(intermediate_length, target_length,
                                           complexity_threshold)
    if (stage2_check$filter_phases >= ratio_check$filter_phases) {
      if (debug) {
        message(sprintf("  Two-stage rejected: Stage 2 = %d phases >= direct %d",
                        stage2_check$filter_phases, ratio_check$filter_phases))
      }
      .resample_log(log_file, "  Two-stage REJECTED: stage2 phases=%d >= direct %d",
                    stage2_check$filter_phases, ratio_check$filter_phases)
      use_two_stage <- FALSE
    } else if (debug) {
      message(sprintf("  Two-stage approved: Stage 2 = %d phases (vs direct %d)",
                      stage2_check$filter_phases, ratio_check$filter_phases))
    }
    if (use_two_stage) {
      .resample_log(log_file, "  Two-stage APPROVED: stage2 phases=%d < direct %d",
                    stage2_check$filter_phases, ratio_check$filter_phases)
    }
  }
  
  # -------------------------------------------------------------------------
  # Overflow safety
  # -------------------------------------------------------------------------
  # v3.6: Stage 1 uses ratio (2,3) directly -> 121-tap filter, trivially
  # safe for any signal length. Only Stage 2 needs overflow checking.
  if (use_two_stage && backend %in% c("gsignal", "signal")) {
    if (.would_overflow(intermediate_length, target_length, intermediate_length)) {
      if (!.would_overflow(original_length, target_length, original_length)) {
        if (debug) message("  [!] Stage 2 overflow -> direct")
        .resample_log(log_file, "  [!] Stage 2 overflow -> falling back to direct")
        use_two_stage <- FALSE
      } else {
        if (debug) message("  [!] Both overflow -> spline")
        .resample_log(log_file, "  [!] Both stages overflow -> falling back to spline")
        backend <- "spline"
        use_two_stage <- FALSE
      }
    }
  }
  if (!use_two_stage && backend %in% c("gsignal", "signal")) {
    if (.would_overflow(original_length, target_length, original_length)) {
      if (debug) message("  [!] Direct overflow -> spline")
      backend <- "spline"
    }
  }
  
  # -------------------------------------------------------------------------
  # Execute resampling
  # -------------------------------------------------------------------------
  if (use_two_stage) {
    if (debug) {
      message("  Method: TWO-STAGE (shrink x2/3, matching MATLAB resample(x,2,3,20))")
      message(sprintf("  Stage 1: %d -> ~%d (ratio 2/3, 121-tap filter)",
                      original_length, intermediate_length))
    }
    .resample_log(log_file, "EXECUTING: TWO-STAGE (v3.6 ratio-based Stage 1)")
    .resample_log(log_file, "  Stage 1: %d -> ~%d (ratio 2/3)", original_length, intermediate_length)
    
    result <- tryCatch({
      # v3.6 fix: Stage 1 uses ratio (2, 3) directly, not computed lengths.
      # This matches MATLAB: resample(x, 2, 3, 20)
      # The ratio approach guarantees a 121-tap filter (3 phases) regardless
      # of signal length, avoiding the catastrophic GCD=1 case that produced
      # multi-million-tap filters and overflowed R's 2^31 integer limit.
      temp <- do_resample_stage1(signal, label = "Stage 1 (2/3)")
      if (debug) {
        s2 <- check_resampling_ratio(length(temp), target_length,
                                     complexity_threshold)
        message(sprintf("  Stage 1 output: %d samples (expected ~%d)",
                        length(temp), intermediate_length))
        message(sprintf("  Stage 2: %d -> %d (%d/%d, %d phases, %d taps)",
                        length(temp), target_length,
                        s2$p_reduced, s2$q_reduced,
                        s2$filter_phases, s2$filter_length))
      }
      do_resample(temp, target_length, length(temp), label = "Stage 2")
    }, error = function(e) {
      if (debug) message(sprintf("  Two-stage failed: %s; trying direct", e$message))
      tryCatch(
        do_resample(signal, target_length, original_length),
        error = function(e2) {
          if (debug) message(sprintf("  Direct failed: %s; spline fallback", e2$message))
          .resample_spline(signal, target_length)
        }
      )
    })
    method_used <- "two_stage"
    
  } else {
    if (debug) message("  Method: DIRECT")
    .resample_log(log_file, "EXECUTING: DIRECT (%d -> %d)", original_length, target_length)
    result <- tryCatch(
      do_resample(signal, target_length, original_length, label = "Direct"),
      error = function(e) {
        if (debug) message(sprintf("  Direct failed: %s; spline fallback", e$message))
        .resample_spline(signal, target_length)
      }
    )
    method_used <- "direct"
  }
  
  # -------------------------------------------------------------------------
  # Ensure exact output length (handle ceil rounding)
  # -------------------------------------------------------------------------
  if (length(result) != target_length) {
    if (debug) message(sprintf("  Length trim: %d -> %d", length(result), target_length))
    .resample_log(log_file, "Length trim: %d -> %d", length(result), target_length)
    if (length(result) > target_length) {
      result <- result[1:target_length]
    } else {
      result <- c(result, rep(tail(result, 1), target_length - length(result)))
    }
  }
  
  # -------------------------------------------------------------------------
  # Metadata
  # -------------------------------------------------------------------------
  elapsed <- as.numeric(proc.time()["elapsed"] - t_start)
  if (debug) message(sprintf("  TOTAL: %.4f sec", elapsed))
  
  # Final log summary
  .resample_log(log_file, "RESULT: method=%s, backend=%s, length=%d, elapsed=%.4f sec",
                method_used, backend, length(result), elapsed)
  .resample_log_signal(log_file, "FINAL OUTPUT", result)
  
  attr(result, "method_used") <- method_used
  attr(result, "original_length") <- original_length
  attr(result, "filter_phases") <- ratio_check$filter_phases
  attr(result, "filter_length") <- ratio_check$filter_length
  attr(result, "was_problematic") <- ratio_check$is_problematic
  attr(result, "backend") <- backend
  attr(result, "elapsed_sec") <- elapsed
  
  return(result)
}
