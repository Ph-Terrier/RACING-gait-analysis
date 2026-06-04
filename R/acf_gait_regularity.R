# =============================================================================
# ACF-Based Gait Regularity Analysis for RACING Project
# =============================================================================
#
# Computes step and stride regularity from trunk accelerometer data using
# the autocorrelation function (ACF) method. Step regularity corresponds to
# the first ACF peak (single step period), while stride regularity corresponds
# to the second peak (full gait cycle).
#
# Scientific References:
#   Primary methodology:
#     Moe-Nilssen, R., & Helbostad, J. L. (2004). Estimation of gait cycle
#     characteristics by trunk accelerometry. Journal of Biomechanics, 37(1), 121-126.
#     https://doi.org/10.1016/S0021-9290(03)00233-1
#
#   Fisher z-transformation for correlation analysis:
#     Fisher, R. A. (1921). On the probable error of a coefficient of correlation
#     deduced from a small sample. Metron, 1, 3-32.
#
#   ACIER study application:
#     Piergiovanni, S., & Terrier, P. (2024). Validity of linear and nonlinear
#     measures of gait variability to characterize aging gait with a single lower
#     back accelerometer. Sensors, 24(23), 7427.
#     https://doi.org/10.3390/s24237427
#
#   Fisher z-transform in gait analysis (recommended applying atanh):
#     Auvinet, B., Berrut, G., Touzard, C., Moutel, L., Collet, N.,
#     Chaleil, D., & Barrey, E. (2002). Reference data for normal subjects
#     obtained with an accelerometric device. Gait & Posture, 16(2), 124-134.
#     https://doi.org/10.1016/S0966-6362(01)00203-X
#
# ACIER Implementation Notes:
# - Input must be the VECTOR NORM of 3D acceleration with gravity removed:
#   norm = sqrt(ax^2 + ay^2 + az^2) - 1
# - The truncated signal (standardized to 250 steps) is expected
# - Fisher z-transform (atanh) is applied for statistical normalization
#
# Author: RACING ORD 2025 Project
# =============================================================================

# -----------------------------------------------------------------------------
# peakn: N-order peak detection (MATLAB equivalent)
# -----------------------------------------------------------------------------
#' Detect local maxima with n-order neighborhood criterion
#'
#' A point is considered a peak only if it is the maximum value within a 
#' window of +/- n samples. This replicates the MATLAB peakn function used
#' in the ACIER study.
#'
#' @param x Numeric vector to search for peaks
#' @param n Order of peak detection (window half-width in samples)
#'
#' @return Data frame with columns:
#'   - index: Position of peak in x (1-indexed)
#'   - value: Value of x at peak position
#'
#' @details
#' The algorithm iterates through valid positions (n+1 to length(x)-n) and
#' checks if each point is the maximum within its neighborhood window.
#' Peaks are returned in order of appearance (by index), NOT by magnitude.
#'
#' @examples
#' x <- c(1, 3, 2, 5, 4, 6, 3, 8, 2)
#' peakn(x, 1)  # Find peaks that dominate +/- 1 neighbor

peakn <- function(x, n) {
  
  # Input validation
  if (!is.numeric(x)) {
    stop("Input x must be a numeric vector")
  }
  if (length(x) < 2 * n + 1) {
    stop(sprintf("Signal too short: need at least %d samples for n=%d", 
                 2 * n + 1, n))
  }
  if (n < 1) {
    stop("Peak order n must be >= 1")
  }
  
  # Initialize output
  peaks <- data.frame(index = integer(), value = numeric())
  
  # Search for peaks in valid range
  for (i in (n + 1):(length(x) - n)) {
    # Extract window centered on position i
    window_start <- i - n
    window_end <- i + n
    window <- x[window_start:window_end]
    
    # Check if center point (position n+1 in window) is the maximum
    max_pos <- which.max(window)
    
    if (max_pos == n + 1) {
      # This point is a local maximum within the n-neighborhood
      peaks <- rbind(peaks, data.frame(index = i, value = x[i]))
    }
  }
  
  return(peaks)
}


# -----------------------------------------------------------------------------
# xcorr_unbiased: Unbiased autocorrelation (MATLAB xcorr equivalent)
# -----------------------------------------------------------------------------
#' Compute unbiased autocorrelation function
#'
#' Replicates MATLAB's xcorr(x, maxlag, 'unbiased') behavior exactly.
#' The 'unbiased' option divides by (n - |k|) at each lag k, compensating
#' for the reduced number of overlapping samples at larger lags.
#'
#' @param x Numeric vector (NOT mean-centered, matching MATLAB behavior)
#' @param max_lag Maximum lag to compute (in samples)
#' @param normalize If TRUE (default), normalize by zero-lag value to get
#'                  correlation coefficients in range [-1, 1]
#'
#' @return Numeric vector of length max_lag + 1, containing autocorrelation
#'         values for lags 0 to max_lag (positive lags only)
#'
#' @details
#' The unbiased estimator is: r(k) = (1/(n-k)) * sum(x[t] * x[t+k])
#' IMPORTANT: MATLAB's xcorr does NOT mean-center the signal. The normalization
#' by r(0) effectively handles centering for correlation purposes.

xcorr_unbiased <- function(x, max_lag, normalize = TRUE) {
  
  # Input validation
  if (!is.numeric(x)) {
    stop("Input x must be a numeric vector")
  }
  n <- length(x)
  if (max_lag >= n) {
    stop(sprintf("max_lag (%d) must be less than signal length (%d)", 
                 max_lag, n))
  }
  

  # NOTE: MATLAB's xcorr does NOT mean-center the signal
  # It computes raw cross-correlation, not correlation coefficients
  # Mean-centering happens implicitly only when normalizing by r(0)
  
  # Compute autocorrelation for lags 0 to max_lag
  r <- numeric(max_lag + 1)
  
  for (k in 0:max_lag) {
    # Sum of products for lag k (NO mean-centering, matching MATLAB)
    sum_prod <- sum(x[1:(n - k)] * x[(k + 1):n])
    # Unbiased normalization: divide by (n - k)
    r[k + 1] <- sum_prod / (n - k)
  }
  
  # Normalize by zero-lag value to get correlation coefficients
  if (normalize) {
    r <- r / r[1]
  }
  
  return(r)
}


# -----------------------------------------------------------------------------
# compute_gait_regularity: Main function
# -----------------------------------------------------------------------------
#' Compute ACF-based gait regularity metrics
#'
#' Analyzes trunk acceleration data using the autocorrelation function to
#' extract step regularity (1st peak) and stride regularity (2nd peak).
#'
#' @param acc_norm Numeric vector: vector magnitude of 3D acceleration with
#'                 gravity removed. Expected input is the truncated signal
#'                 from the ACIER preprocessing pipeline.
#'                 IMPORTANT: Must be computed as sqrt(ax^2 + ay^2 + az^2) - 1
#' @param sampling_freq Sampling frequency in Hz (default: 256)
#' @param max_lag_seconds Maximum lag for ACF computation in seconds (default: 2.0)
#' @param peak_order Order for peak detection in samples (default: 60).
#'                   At 256 Hz, 60 samples ~= 234 ms, preventing false peaks
#'                   within the same gait event.
#' @param apply_fisher If TRUE (default), apply Fisher z-transform (atanh) to
#'                     regularity values for statistical normalization
#'
#' @return List containing:
#'   \item{step_regularity}{First ACF peak value (Fisher-transformed if requested)}
#'   \item{stride_regularity}{Second ACF peak value (Fisher-transformed if requested)}
#'   \item{step_regularity_raw}{First ACF peak value (untransformed correlation)}
#'   \item{stride_regularity_raw}{Second ACF peak value (untransformed correlation)}
#'   \item{step_lag_samples}{Lag of first peak in samples}
#'   \item{stride_lag_samples}{Lag of second peak in samples}
#'   \item{step_period_seconds}{Estimated step period in seconds}
#'   \item{stride_period_seconds}{Estimated stride period in seconds}
#'   \item{acf_values}{Full ACF sequence (for plotting)}
#'   \item{peaks}{Data frame of all detected peaks}
#'   \item{parameters}{List of input parameters used}
#'
#' @details
#' The algorithm follows the ACIER/MATLAB implementation:
#' 1. Compute unbiased autocorrelation up to max_lag
#' 2. Normalize by zero-lag value (correlation coefficients)
#' 3. Detect peaks using n-order criterion (peakn)
#' 4. Extract 1st peak (step) and 2nd peak (stride) values
#' 5. Apply Fisher z-transform for statistical analysis
#'
#' @references
#' Moe-Nilssen R, Helbostad JL (2004). Estimation of gait cycle characteristics
#' by trunk accelerometry. J Biomech, 37(1):121-6.
#'
#' Auvinet B, Berrut G, Touzard C, et al. (2002). Reference data for normal
#' subjects obtained with an accelerometric device. Gait Posture, 16(2):124-34.
#'
#' @examples
#' # Assuming acc_norm is the preprocessed vector magnitude signal:
#' # result <- compute_gait_regularity(acc_norm)
#' # print(result$step_regularity)
#' # print(result$stride_regularity)

compute_gait_regularity <- function(acc_norm,
                                     sampling_freq = 256,
                                     max_lag_seconds = 2.0,
                                     peak_order = 60,
                                     apply_fisher = TRUE) {
  
 
  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if (!is.numeric(acc_norm)) {
    stop("acc_norm must be a numeric vector")
  }
  
  if (any(is.na(acc_norm))) {
    warning("Input contains NA values - removing them")
    acc_norm <- acc_norm[!is.na(acc_norm)]
  }
  
  n_samples <- length(acc_norm)
  max_lag_samples <- floor(max_lag_seconds * sampling_freq)
  
  # Check signal length
  min_required <- max_lag_samples + peak_order + 1
  if (n_samples < min_required) {
    stop(sprintf(
      "Signal too short: %d samples, need at least %d for max_lag=%.1fs and peak_order=%d",
      n_samples, min_required, max_lag_seconds, peak_order
    ))
  }
  
  # -------------------------------------------------------------------------
  # Step 1: Compute unbiased autocorrelation
  # -------------------------------------------------------------------------
  acf_values <- xcorr_unbiased(acc_norm, max_lag_samples, normalize = TRUE)
  
  # -------------------------------------------------------------------------
  # Step 2: Detect peaks using n-order criterion
  # -------------------------------------------------------------------------
  peaks <- peakn(acf_values, peak_order)
  
  # Check if we found at least 2 peaks
  if (nrow(peaks) < 2) {
    warning(sprintf(
      "Only %d peak(s) detected. Need at least 2 for step and stride regularity. %s",
      nrow(peaks),
      "Consider reducing peak_order or checking signal quality."
    ))
    
    # Return NA values if insufficient peaks
    return(list(
      step_regularity = NA_real_,
      stride_regularity = NA_real_,
      step_regularity_raw = NA_real_,
      stride_regularity_raw = NA_real_,
      step_lag_samples = NA_integer_,
      stride_lag_samples = NA_integer_,
      step_period_seconds = NA_real_,
      stride_period_seconds = NA_real_,
      acf_values = acf_values,
      peaks = peaks,
      parameters = list(
        sampling_freq = sampling_freq,
        max_lag_seconds = max_lag_seconds,
        max_lag_samples = max_lag_samples,
        peak_order = peak_order,
        apply_fisher = apply_fisher,
        n_samples = n_samples
      )
    ))
  }
  
  # -------------------------------------------------------------------------
  # Step 3: Extract step (1st peak) and stride (2nd peak) values
  # -------------------------------------------------------------------------
  # Peaks are returned in order of appearance (by lag)
  step_peak <- peaks[1, ]
  stride_peak <- peaks[2, ]
  
  # Raw correlation values
  step_reg_raw <- step_peak$value
  stride_reg_raw <- stride_peak$value
  
  # Lag positions (convert from 1-indexed R to 0-indexed lag)
  step_lag <- step_peak$index - 1
  stride_lag <- stride_peak$index - 1
  
  # -------------------------------------------------------------------------
  # Step 4: Apply Fisher z-transform if requested
  # -------------------------------------------------------------------------
  if (apply_fisher) {
    # Clamp values to valid range for atanh (-1, 1) exclusive
    # This handles potential numerical issues
    step_reg_clamped <- max(-0.9999, min(0.9999, step_reg_raw))
    stride_reg_clamped <- max(-0.9999, min(0.9999, stride_reg_raw))
    
    step_regularity <- atanh(step_reg_clamped)
    stride_regularity <- atanh(stride_reg_clamped)
  } else {
    step_regularity <- step_reg_raw
    stride_regularity <- stride_reg_raw
  }
  
  # -------------------------------------------------------------------------
  # Step 5: Compute period estimates
  # -------------------------------------------------------------------------
  step_period <- step_lag / sampling_freq
  stride_period <- stride_lag / sampling_freq
  
  # -------------------------------------------------------------------------
  # Return results
  # -------------------------------------------------------------------------
  return(list(
    step_regularity = step_regularity,
    stride_regularity = stride_regularity,
    step_regularity_raw = step_reg_raw,
    stride_regularity_raw = stride_reg_raw,
    step_lag_samples = step_lag,
    stride_lag_samples = stride_lag,
    step_period_seconds = step_period,
    stride_period_seconds = stride_period,
    acf_values = acf_values,
    peaks = peaks,
    parameters = list(
      sampling_freq = sampling_freq,
      max_lag_seconds = max_lag_seconds,
      max_lag_samples = max_lag_samples,
      peak_order = peak_order,
      apply_fisher = apply_fisher,
      n_samples = n_samples
    )
  ))
}


# -----------------------------------------------------------------------------
# plot_gait_regularity: Visualization function
# -----------------------------------------------------------------------------
#' Plot ACF with detected step and stride peaks
#'
#' Creates a visualization of the autocorrelation function with markers
#' indicating the detected step period (1st peak) and stride period (2nd peak).
#'
#' @param result Output from compute_gait_regularity()
#' @param show_all_peaks If TRUE, mark all detected peaks (default: FALSE)
#' @param title Optional custom title for the plot
#' @param time_axis If TRUE (default), show x-axis in seconds; if FALSE, in samples
#'
#' @return Invisibly returns the result object (for chaining)
#'
#' @examples
#' # result <- compute_gait_regularity(acc_norm)
#' # plot_gait_regularity(result)

plot_gait_regularity <- function(result,
                                  show_all_peaks = FALSE,
                                  title = NULL,
                                  time_axis = TRUE) {
  
  # Check input
 if (!is.list(result) || !"acf_values" %in% names(result)) {
    stop("Input must be output from compute_gait_regularity()")
  }
  
  # Save and restore graphics parameters
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Set up plotting parameters for RStudio compatibility
  par(mar = c(4, 4, 3, 1),   # Margins: bottom, left, top, right
      mgp = c(2.5, 0.8, 0),  # Axis label positions
      cex = 0.9)             # Text scaling
  
  # Extract data
  acf_vals <- result$acf_values
  peaks <- result$peaks
  params <- result$parameters
  
  # Create x-axis values
  n_lags <- length(acf_vals)
  if (time_axis) {
    x_vals <- (0:(n_lags - 1)) / params$sampling_freq
    x_label <- "Lag (seconds)"
  } else {
    x_vals <- 0:(n_lags - 1)
    x_label <- "Lag (samples)"
  }
  
  # Set title
  if (is.null(title)) {
    title <- "ACF-Based Gait Regularity Analysis"
  }
  
  # Create plot
  plot(x_vals, acf_vals,
       type = "l",
       col = "darkblue",
       lwd = 1.5,
       xlab = x_label,
       ylab = "Autocorrelation",
       main = title,
       ylim = c(min(acf_vals, -0.2), 1.0))
  
  # Add horizontal reference line at 0
  abline(h = 0, col = "gray50", lty = 2)
  
  # Mark all peaks if requested
  if (show_all_peaks && nrow(peaks) > 0) {
    if (time_axis) {
      peak_x <- (peaks$index - 1) / params$sampling_freq
    } else {
      peak_x <- peaks$index - 1
    }
    points(peak_x, peaks$value,
           pch = 1, col = "gray40", cex = 1.2)
  }
  
  # Mark step peak (1st peak) - green
  if (!is.na(result$step_lag_samples)) {
    if (time_axis) {
      step_x <- result$step_period_seconds
    } else {
      step_x <- result$step_lag_samples
    }
    points(step_x, result$step_regularity_raw,
           pch = 19, col = "forestgreen", cex = 2)
    
    # Add label
    text(step_x, result$step_regularity_raw,
         labels = sprintf("Step\n(%.3f)", result$step_regularity_raw),
         pos = 4, col = "forestgreen", cex = 0.8)
  }
  
  # Mark stride peak (2nd peak) - red
  if (!is.na(result$stride_lag_samples)) {
    if (time_axis) {
      stride_x <- result$stride_period_seconds
    } else {
      stride_x <- result$stride_lag_samples
    }
    points(stride_x, result$stride_regularity_raw,
           pch = 19, col = "firebrick", cex = 2)
    
    # Add label
    text(stride_x, result$stride_regularity_raw,
         labels = sprintf("Stride\n(%.3f)", result$stride_regularity_raw),
         pos = 4, col = "firebrick", cex = 0.8)
  }
  
  # Add legend with regularity values
  if (!is.na(result$step_regularity)) {
    legend_text <- c(
      sprintf("Step: %.3f (z=%.3f) | T=%.3f s", 
              result$step_regularity_raw, result$step_regularity, 
              result$step_period_seconds),
      sprintf("Stride: %.3f (z=%.3f) | T=%.3f s", 
              result$stride_regularity_raw, result$stride_regularity,
              result$stride_period_seconds)
    )
    
    legend("topright",
           legend = legend_text,
           col = c("forestgreen", "firebrick"),
           pch = c(19, 19),
           cex = 0.75,
           bg = "white")
  }
  
  # Return result invisibly for chaining
  invisible(result)
}


# -----------------------------------------------------------------------------
# print method for gait regularity results
# -----------------------------------------------------------------------------
#' Print summary of gait regularity results
#'
#' @param result Output from compute_gait_regularity()

print_gait_regularity <- function(result) {
  
  cat("\n========================================\n")
  cat("  ACF Gait Regularity Results\n")
  cat("========================================\n\n")
  
  if (is.na(result$step_regularity)) {
    cat("  WARNING: Insufficient peaks detected\n\n")
  }
  
  cat("Step Regularity:\n")
  cat(sprintf("  Raw correlation:    %.4f\n", result$step_regularity_raw))
  cat(sprintf("  Fisher z-transform: %.4f\n", result$step_regularity))
  cat(sprintf("  Period:             %.3f s (%d samples)\n", 
              result$step_period_seconds, result$step_lag_samples))
  
  cat("\nStride Regularity:\n")
  cat(sprintf("  Raw correlation:    %.4f\n", result$stride_regularity_raw))
  cat(sprintf("  Fisher z-transform: %.4f\n", result$stride_regularity))
  cat(sprintf("  Period:             %.3f s (%d samples)\n", 
              result$stride_period_seconds, result$stride_lag_samples))
  
  cat("\nDetected Peaks:\n")
  cat(sprintf("  Total peaks found: %d\n", nrow(result$peaks)))
  
  cat("\nParameters Used:\n")
  cat(sprintf("  Sampling frequency: %d Hz\n", result$parameters$sampling_freq))
  cat(sprintf("  Max lag:            %.2f s (%d samples)\n", 
              result$parameters$max_lag_seconds, result$parameters$max_lag_samples))
  cat(sprintf("  Peak order:         %d samples\n", result$parameters$peak_order))
  cat(sprintf("  Signal length:      %d samples\n", result$parameters$n_samples))
  
  cat("\n========================================\n")
  
  invisible(result)
}


