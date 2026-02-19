#' =============================================================================
#' Divergence Exponent Fitting from Rosenstein Divergence Curves
#' =============================================================================
#' 
#' This script computes short-term (LDS) and long-term (ACI) divergence exponents
#' from divergence curves produced by the rosenstein_divergence() function.
#' 
#' ## TIME NORMALIZATION
#' 
#' The divergence curve x-axis is normalized by STRIDE duration (not step):
#' - 1 STEP = 75 samples (in ACIER standardized signals at 256 Hz)
#' - 1 STRIDE = 2 STEPS = 150 samples (~586 ms)
#' - Time axis: t(i) = (i-1) / samples_per_stride, where i is the sample index
#' 
#' ## DIVERGENCE EXPONENT DEFINITIONS
#' 
#' **Short-term Divergence Exponent (LDS - Local Dynamic Stability)**
#' - Measures local instability / sensitivity to perturbations
#' - Computed over the FIRST STEP (0 to 0.5 strides = samples 1 to 75)
#' - Higher values = less stable gait, greater fall risk
#' - Primarily reflects mediolateral (ML) control in frontal plane
#' - Consistent marker of fall risk across multiple studies
#' 
#' **Long-term Divergence Exponent (ACI - Attractor Complexity Index)**
#' - Originally designed as stability measure (2000-2013), reinterpreted as COMPLEXITY (2016+)
#' - Computed over 5-12 strides (samples 750 to 1800)
#' - Correlates strongly with DFA scaling exponents (r ≈ 0.7-0.8)
#' - Reflects stride-to-stride fluctuation patterns and gait automaticity
#' - Lower values during metronome walking indicate reduced automaticity
#' - Not a local stability measure (conceptually distinct from LDS)
#' 
#' ## SCIENTIFIC REFERENCES
#'
#' Rosenstein algorithm (divergence computation):
#'   Rosenstein, M. T., Collins, J. J., & De Luca, C. J. (1993).
#'   A practical method for calculating largest Lyapunov exponents from small data sets.
#'   Physica D: Nonlinear Phenomena, 65(1-2), 117-134.
#'   https://doi.org/10.1016/0167-2789(93)90009-P
#'
#' Short-term divergence exponent (LDS) - gait stability:
#'   Dingwell, J. B., & Cusumano, J. P. (2000).
#'   Nonlinear time series analysis of normal and pathological human walking.
#'   Chaos, 10(4), 848-863.
#'   https://doi.org/10.1063/1.1324008
#'
#'   Terrier, P., & Reynard, F. (2018).
#'   Maximum Lyapunov exponent revisited: Long-term attractor divergence of gait 
#'   dynamics is highly sensitive to the noise structure of stride intervals.
#'   Gait & Posture, 66, 236-241.
#'   https://doi.org/10.1016/j.gaitpost.2018.08.010
#'
#' ACIER study - empirical validation and age effects:
#'   Piergiovanni, S., & Terrier, P. (2024).
#'   Effects of metronome walking on long-term attractor divergence and correlation
#'   structure of gait: A validation study in older people.
#'   Scientific Reports, 14, 15784.
#'   https://doi.org/10.1038/s41598-024-65662-5
#'
#'   Piergiovanni, S., & Terrier, P. (2024).
#'   Validity of linear and nonlinear measures of gait variability to characterize
#'   aging gait with a single lower back accelerometer.
#'   Sensors, 24(23), 7427.
#'   https://doi.org/10.3390/s24237427
#'
#' ## IMPLEMENTATION NOTES
#'
#' This R implementation replicates the MATLAB code from the ACIER study
#' (Main_results_new.m), including:
#' - Identical fitting ranges (samples 1:75 for LDS, 750:1800 for ACI)
#' - Same time normalization (by stride duration)
#' - Linear regression via ordinary least squares
#' - Stride range optimization (5-12 strides) from systematic sensitivity analysis
#'
#' Author: Philippe Terrier, RACING Project (ORD 2025), HE-Arc 
#' 
#' =============================================================================
#' Fit Short-term Divergence Exponent (LDS)
#' 
#' Computes the Local Dynamic Stability index by fitting a linear regression
#' to the initial portion of the divergence curve (first step = 0-0.5 strides).
#' 
#' @param divergence Numeric vector - divergence curve from rosenstein_divergence()
#' @param samples_per_step Integer - number of samples per step (default: 75 for ACIER)
#' @param num_steps Integer - number of steps for LDS fitting (default: 1)
#' @param return_fit Logical - if TRUE, return full lm object (default: FALSE)
#' 
#' @return If return_fit=FALSE: numeric value (LDS slope in 1/stride units)
#'         If return_fit=TRUE: list with slope, intercept, r_squared, fit object
#' 
#' @details
#' Time is normalized by STRIDE (2 steps), so:
#' - 1 step = samples_per_step samples = 0.5 strides
#' - Default fitting range: indices 1 to 75 = time 0 to ~0.493 strides
#' 
#' The slope represents the average exponential divergence rate:
#' d(t) = d0 * exp(lambda * t), so log(d(t)) ≈ log(d0) + lambda * t
#' 
#' @examples
#' result <- rosenstein_divergence(gait_signal, m=5, tau=10, max_iter=1800)
#' lds <- fit_short_term_divergence(result$divergence)
fit_short_term_divergence <- function(divergence, 
                                       samples_per_step = 75,
                                       num_steps = 1,
                                       return_fit = FALSE) {
  
  # Validate inputs
  if (!is.numeric(divergence)) {
    stop("divergence must be a numeric vector")
  }
  if (samples_per_step <= 0 || !is.numeric(samples_per_step)) {
    stop("samples_per_step must be a positive number")
  }
  if (num_steps <= 0) {
    stop("num_steps must be positive")
  }
  
  # Calculate fitting range (MATLAB: 1:numsamp*lds_length)
  # In R, indices are 1-based like MATLAB
  end_idx <- samples_per_step * num_steps
  
  # Check we have enough data
  if (end_idx > length(divergence)) {
    stop(sprintf("Divergence curve too short (%d samples). Need at least %d for %d step(s).",
                 length(divergence), end_idx, num_steps))
  }
  
  # Extract fitting range
  indices <- 1:end_idx
  div_segment <- divergence[indices]
  
  # Check for valid data
  valid_mask <- !is.na(div_segment) & is.finite(div_segment)
  if (sum(valid_mask) < 2) {
    warning("Insufficient valid data points for LDS fitting")
    if (return_fit) {
      return(list(slope = NA, intercept = NA, r_squared = NA, fit = NULL,
                  range_samples = c(1, end_idx), range_strides = c(0, num_steps/2)))
    }
    return(NA)
  }
  
  # Create time axis normalized by STRIDE (2 steps = samples_per_step * 2 samples)
  # MATLAB: t = 0:1/(numsamp*2):div_length/(numsamp*2)
  # For index i (1-based), time = (i-1) / (samples_per_step * 2)
  samples_per_stride <- samples_per_step * 2
  time_strides <- (indices - 1) / samples_per_stride
  
  # Linear fit using only valid points
  fit_data <- data.frame(time = time_strides[valid_mask], 
                         div = div_segment[valid_mask])
  fit <- lm(div ~ time, data = fit_data)
  
  slope <- coef(fit)[2]  # Lambda_s (short-term)
  
  if (return_fit) {
    return(list(
      slope = as.numeric(slope),
      intercept = as.numeric(coef(fit)[1]),
      r_squared = summary(fit)$r.squared,
      fit = fit,
      range_samples = c(1, end_idx),
      range_strides = c(0, (end_idx - 1) / samples_per_stride),
      n_points = sum(valid_mask)
    ))
  }
  
  as.numeric(slope)
}


#' Fit Long-term Divergence Exponent (ACI)
#' 
#' Computes the Attractor Complexity Index by fitting a linear regression
#' to the later portion of the divergence curve (typically 5-12 strides).
#' 
#' @param divergence Numeric vector - divergence curve from rosenstein_divergence()
#' @param samples_per_step Integer - number of samples per step (default: 75 for ACIER)
#' @param stride_range Numeric vector of length 2 - fitting range in strides (default: c(5, 12))
#' @param return_fit Logical - if TRUE, return full lm object (default: FALSE)
#' 
#' @return If return_fit=FALSE: numeric value (ACI slope in 1/stride units)
#'         If return_fit=TRUE: list with slope, intercept, r_squared, fit object
#' 
#' @details
#' ACIER MATLAB implementation uses strides 5-12:
#' - Start index: 5 * 75 * 2 = 750
#' - End index: 12 * 75 * 2 = 1800
#' 
#' The ACI reflects long-term attractor dynamics:
#' - Lower ACI = faster divergence saturation = reduced complexity
#' - Correlates with DFA alpha (fractal index of stride intervals)
#' - NOT a stability measure (despite being derived from divergence)
#' 
#' @examples
#' result <- rosenstein_divergence(gait_signal, m=5, tau=10, max_iter=1800)
#' aci <- fit_long_term_divergence(result$divergence)
fit_long_term_divergence <- function(divergence,
                                      samples_per_step = 75,
                                      stride_range = c(5, 12),
                                      return_fit = FALSE) {
  
  # Validate inputs
  if (!is.numeric(divergence)) {
    stop("divergence must be a numeric vector")
  }
  if (samples_per_step <= 0) {
    stop("samples_per_step must be positive")
  }
  if (length(stride_range) != 2 || stride_range[1] >= stride_range[2]) {
    stop("stride_range must be c(start, end) with start < end")
  }
  if (stride_range[1] < 0) {
    stop("stride_range values must be non-negative")
  }
  
  # Calculate fitting range (MATLAB: range_low*numsamp*2:range_high*numsamp*2)
  samples_per_stride <- samples_per_step * 2
  start_idx <- stride_range[1] * samples_per_stride
  end_idx <- stride_range[2] * samples_per_stride
  
  # MATLAB uses 1-based indexing, but the time starts at 0
  # Index 750 corresponds to time = 749/150 ≈ 4.99 strides
  # We need to match MATLAB exactly: indices from start_idx to end_idx
  # In R, this is the same as MATLAB (1-based)
  
  # Check we have enough data
  if (end_idx > length(divergence)) {
    stop(sprintf("Divergence curve too short (%d samples). Need at least %d for %.1f-%.1f stride range.",
                 length(divergence), end_idx, stride_range[1], stride_range[2]))
  }
  
  # Extract fitting range
  indices <- start_idx:end_idx
  div_segment <- divergence[indices]
  
  # Check for valid data
  valid_mask <- !is.na(div_segment) & is.finite(div_segment)
  if (sum(valid_mask) < 2) {
    warning("Insufficient valid data points for ACI fitting")
    if (return_fit) {
      return(list(slope = NA, intercept = NA, r_squared = NA, fit = NULL,
                  range_samples = c(start_idx, end_idx), 
                  range_strides = stride_range))
    }
    return(NA)
  }
  
  # Create time axis normalized by STRIDE
  # For index i (1-based), time = (i-1) / samples_per_stride
  time_strides <- (indices - 1) / samples_per_stride
  
  # Linear fit using only valid points
  fit_data <- data.frame(time = time_strides[valid_mask], 
                         div = div_segment[valid_mask])
  fit <- lm(div ~ time, data = fit_data)
  
  slope <- coef(fit)[2]  # Lambda_l (long-term) = ACI
  
  if (return_fit) {
    return(list(
      slope = as.numeric(slope),
      intercept = as.numeric(coef(fit)[1]),
      r_squared = summary(fit)$r.squared,
      fit = fit,
      range_samples = c(start_idx, end_idx),
      range_strides = c((start_idx - 1) / samples_per_stride, 
                        (end_idx - 1) / samples_per_stride),
      n_points = sum(valid_mask)
    ))
  }
  
  as.numeric(slope)
}


#' Fit Both Divergence Exponents (LDS and ACI)
#' 
#' Convenience function to compute both short-term (LDS) and long-term (ACI)
#' divergence exponents from a single divergence curve.
#' 
#' @param divergence Numeric vector - divergence curve from rosenstein_divergence()
#' @param samples_per_step Integer - number of samples per step (default: 75)
#' @param lds_steps Integer - number of steps for LDS fitting (default: 1)
#' @param aci_stride_range Numeric vector - stride range for ACI (default: c(5, 12))
#' @param return_fits Logical - if TRUE, include full lm objects (default: FALSE)
#' 
#' @return List containing:
#'   - LDS: short-term divergence exponent (1/stride)
#'   - ACI: long-term divergence exponent / Attractor Complexity Index (1/stride)
#'   - params: analysis parameters used
#'   - fits: (optional) detailed fit information for both
#' 
#' @examples
#' # Standard ACIER analysis
#' result <- rosenstein_divergence(signal, m=5, tau=10, max_iter=1800)
#' exponents <- fit_divergence_exponents(result$divergence)
#' 
#' # Custom ranges
#' exponents <- fit_divergence_exponents(result$divergence,
#'                                        lds_steps = 1,
#'                                        aci_stride_range = c(4, 10))
fit_divergence_exponents <- function(divergence,
                                      samples_per_step = 75,
                                      lds_steps = 1,
                                      aci_stride_range = c(5, 12),
                                      return_fits = FALSE) {
  
  # Compute LDS
  lds_result <- fit_short_term_divergence(divergence, 
                                           samples_per_step = samples_per_step,
                                           num_steps = lds_steps,
                                           return_fit = return_fits)
  
  # Compute ACI
  aci_result <- fit_long_term_divergence(divergence,
                                          samples_per_step = samples_per_step,
                                          stride_range = aci_stride_range,
                                          return_fit = return_fits)
  
  # Build output
  samples_per_stride <- samples_per_step * 2
  
  output <- list(
    LDS = if (return_fits) lds_result$slope else lds_result,
    ACI = if (return_fits) aci_result$slope else aci_result,
    params = list(
      samples_per_step = samples_per_step,
      samples_per_stride = samples_per_stride,
      lds_steps = lds_steps,
      lds_range_samples = c(1, samples_per_step * lds_steps),
      lds_range_strides = c(0, lds_steps / 2),
      aci_stride_range = aci_stride_range,
      aci_range_samples = c(aci_stride_range[1] * samples_per_stride,
                            aci_stride_range[2] * samples_per_stride)
    )
  )
  
  if (return_fits) {
    output$fits <- list(
      LDS = lds_result,
      ACI = aci_result
    )
  }
  
  output
}


#' Plot Divergence Curve with Fitted Exponents
#' 
#' Visualizes the divergence curve with LDS and ACI fitting regions highlighted.
#' 
#' @param divergence Numeric vector - divergence curve
#' @param samples_per_step Integer - samples per step (default: 75)
#' @param lds_steps Integer - steps for LDS (default: 1)
#' @param aci_stride_range Numeric vector - strides for ACI (default: c(5, 12))
#' @param show_legend Logical - display legend (default: TRUE)
#' @param main Character - plot title (default: "Divergence Curve Analysis")
#' @param ... Additional arguments passed to plot()
#' 
#' @return Invisibly returns the fitted exponents list
#' 
#' @examples
#' plot_divergence_fits(divergence_curve)
plot_divergence_fits <- function(divergence,
                                  samples_per_step = 75,
                                  lds_steps = 1,
                                  aci_stride_range = c(5, 12),
                                  show_legend = TRUE,
                                  main = "Divergence Curve Analysis",
                                  ...) {
  
  # Save/restore graphics parameters (critical for RStudio)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Use RStudio-compatible margins
  par(mar = c(2, 2, 1.5, 0.5), 
      mgp = c(1.8, 0.5, 0),
      cex = 0.7)
  
  # Get fitted values
  fits <- fit_divergence_exponents(divergence, 
                                    samples_per_step = samples_per_step,
                                    lds_steps = lds_steps,
                                    aci_stride_range = aci_stride_range,
                                    return_fits = TRUE)
  
  # Create time axis in strides
  samples_per_stride <- samples_per_step * 2
  n <- length(divergence)
  time_strides <- (seq_len(n) - 1) / samples_per_stride
  
  # Remove NA values for plotting
  valid <- !is.na(divergence) & is.finite(divergence)
  
  # Plot main curve
  plot(time_strides[valid], divergence[valid], type = "l", 
       xlab = "Time (strides)", ylab = expression("<ln d(i)>"),
       main = main, col = "gray40", lwd = 1.5, ...)
  
  # Highlight LDS region and fit
  lds_idx <- 1:(samples_per_step * lds_steps)
  lds_time <- (lds_idx - 1) / samples_per_stride
  lds_valid <- lds_idx[valid[lds_idx]]
  
  if (length(lds_valid) > 0) {
    lines(time_strides[lds_valid], divergence[lds_valid], col = "blue", lwd = 2)
    
    # Add LDS fit line
    if (!is.na(fits$LDS)) {
      lds_fit <- fits$fits$LDS
      lds_pred <- lds_fit$intercept + lds_fit$slope * lds_time
      lines(lds_time, lds_pred, col = "blue", lwd = 2, lty = 2)
    }
  }
  
  # Highlight ACI region and fit
  aci_start <- aci_stride_range[1] * samples_per_stride
  aci_end <- min(aci_stride_range[2] * samples_per_stride, n)
  aci_idx <- aci_start:aci_end
  aci_time <- (aci_idx - 1) / samples_per_stride
  aci_valid <- aci_idx[valid[aci_idx]]
  
  if (length(aci_valid) > 0) {
    lines(time_strides[aci_valid], divergence[aci_valid], col = "red", lwd = 2)
    
    # Add ACI fit line
    if (!is.na(fits$ACI)) {
      aci_fit <- fits$fits$ACI
      aci_pred <- aci_fit$intercept + aci_fit$slope * aci_time
      lines(aci_time, aci_pred, col = "red", lwd = 2, lty = 2)
    }
  }
  
  # Add vertical lines for boundaries
  abline(v = lds_steps / 2, col = "blue", lty = 3)
  abline(v = aci_stride_range[1], col = "red", lty = 3)
  abline(v = aci_stride_range[2], col = "red", lty = 3)
  
  # Legend
  if (show_legend) {
    legend_text <- c(
      sprintf("LDS = %.3f (0-%.1f strides)", fits$LDS, lds_steps/2),
      sprintf("ACI = %.3f (%.0f-%.0f strides)", fits$ACI, 
              aci_stride_range[1], aci_stride_range[2])
    )
    legend("bottomright", legend = legend_text, 
           col = c("blue", "red"), lwd = 2, lty = c(2, 2),
           bty = "n", cex = 0.9)
  }
  
  invisible(fits)
}

 

#' @examples
#' # Standard ACIER analysis pipeline:
#' # 1. Compute divergence curve (from rosenstein_divergence.R)
#' result <- rosenstein_divergence(
#'   gait_signal_tilt_corrected,  # After preprocessing
#'   m = 5,                        # Embedding dimension (standard)
#'   tau = 15,                     # Time delay from AMI analysis
#'   mean_period = 75,             # Temporal exclusion = 1 step
#'   max_iter = 1800               # Must cover ACI range (12 strides)
#' )
#' 
#' # 2. Fit both exponents with default ACIER parameters
#' exponents <- fit_divergence_exponents(result$divergence)
#' cat("LDS (stability):", exponents$LDS, "per stride\n")
#' cat("ACI (complexity):", exponents$ACI, "per stride\n")
#' 
#' # 3. Detailed analysis with fit quality metrics
#' detailed <- fit_divergence_exponents(result$divergence, return_fits = TRUE)
#' cat("LDS R² (should be > 0.95):", detailed$fits$LDS$r_squared, "\n")
#' cat("ACI R² (should be > 0.90):", detailed$fits$ACI$r_squared, "\n")
#' 
#' # 4. Custom ranges for sensitivity analysis
#' custom <- fit_divergence_exponents(
#'   result$divergence,
#'   lds_steps = 1,              # Can test 1, 2, or 3 steps
#'   aci_stride_range = c(4, 10) # Can vary range
#' )
#' 
#' @seealso
#' \code{\link{rosenstein_divergence}} for computing divergence curves
#' \code{\link{plot_divergence_fits}} for visualization