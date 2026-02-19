# =============================================================================
# Step Frequency Detection - FFT-Based Method
# =============================================================================
#
# RACING Project (ORD 2025) - Version 2.1
# Author: Philippe Terrier, HE-Arc Engineering
# Date: February 2026
#
# Purpose: Estimate step frequency from vertical acceleration using Fast
# Fourier Transform (FFT) power spectral density analysis.
#
# Method: Standard FFT-based frequency detection with physiologically
# constrained search range (1.4-3.0 Hz for normal human walking).
#
# Reference: This implements standard spectral analysis methods as described in:
#   Welch, P. (1967). The use of fast Fourier transform for the estimation
#   of power spectra. IEEE Transactions on Audio and Electroacoustics, 15(2), 70-73.
#
# Note: The frequency range (1.4-3.0 Hz) corresponds to 84-180 steps/min,
# covering slow to very fast walking cadences typically observed in gait studies.
#
# =============================================================================
# =============================================================================
# MAIN FUNCTION: FFT-Based Step Frequency Detection
# =============================================================================

#' Compute Step Frequency Using Fast Fourier Transform
#'
#' Estimates step frequency from vertical acceleration using FFT power spectral
#' density analysis. This method matches the ACIER original implementation.
#'
#' @param acceleration_vertical Numeric vector of vertical acceleration (g units)
#' @param sampling_freq Sampling frequency in Hz (default: 256)
#' @param freq_range Physiological frequency range c(min, max) in Hz.
#'   Default: c(1.4, 3.0) covers normal walking (84-180 steps/min).
#' @param detrend Logical. Remove mean from signal (default: TRUE)
#' @param zero_pad Logical. Apply zero-padding for better frequency resolution (default: FALSE)
#' @param return_spectrum Logical. Return full PSD for plotting (default: FALSE)
#'
#' @return List containing:
#'   \item{step_freq_hz}{Step frequency in Hz}
#'   \item{step_freq_bpm}{Step frequency in beats per minute}
#'   \item{peak_frequency}{Frequency of maximum PSD peak}
#'   \item{peak_power}{Power spectral density at peak}
#'   \item{psd}{Power spectral density (if return_spectrum=TRUE)}
#'   \item{frequencies}{Frequency vector (if return_spectrum=TRUE)}
#'
#' @details
#' Algorithm steps:
#' 1. Preprocessing: detrending (remove mean) and optional bandpass filtering
#' 2. FFT computation with optional zero-padding
#' 3. Power spectral density calculation (magnitude squared)
#' 4. Peak detection within physiological frequency range (1.4-3.0 Hz)
#'
#' Frequency range rationale:
#' - 1.4 Hz = 84 steps/min (slow walking, elderly)
#' - 3.0 Hz = 180 steps/min (fast walking)
#'
#' Minimum signal length: 4 seconds (1024 samples at 256 Hz)
#'
#' @examples
#' # Simulate walking signal at 1.8 Hz (108 steps/min)
#' t <- seq(0, 10, 1/256)  # 10 seconds at 256 Hz
#' signal <- sin(2 * pi * 1.8 * t) + 0.2 * rnorm(length(t))
#' 
#' # Compute step frequency
#' result <- step_frequency_fft(signal, sampling_freq = 256)
#' cat("Detected frequency:", round(result$step_freq_hz, 2), "Hz\n")
#' cat("Cadence:", round(result$step_freq_bpm, 1), "steps/min\n")
#'
#' @export
step_frequency_fft <- function(
    acceleration_vertical,
    sampling_freq = 256,
    freq_range = c(1.4, 3.0),
    detrend = TRUE,
    zero_pad = FALSE,
    return_spectrum = FALSE
) {
  
  # --- Input validation ---
  if (!is.numeric(acceleration_vertical) || length(acceleration_vertical) < 100) {
    stop("acceleration_vertical must be a numeric vector with at least 100 samples")
  }
  
  if (length(acceleration_vertical) < 4 * sampling_freq) {
    warning(sprintf(
      "Signal length is %.1f seconds. For best accuracy, use signals >= 4 seconds.",
      length(acceleration_vertical) / sampling_freq
    ))
  }
  
  # --- Preprocessing ---
  signal <- acceleration_vertical
  
  # Remove DC component (detrending)
  if (detrend) {
    signal <- signal - mean(signal)
  }
  
 
  
  # --- FFT computation ---
  N <- length(signal)
  
  # Zero padding for better frequency resolution
  if (zero_pad) {
    N_padded <- 2^ceiling(log2(N))
    signal_padded <- c(signal, rep(0, N_padded - N))
  } else {
    signal_padded <- signal
    N_padded <- N
  }
  
  # Compute FFT
  fft_result <- fft(signal_padded)
  
  # Take only positive frequencies (first half + DC)
  N_half <- floor(N_padded / 2) + 1
  fft_positive <- fft_result[1:N_half]
  
  # --- Power Spectral Density ---
  # One-sided power spectral density: magnitude squared, normalized by signal length
  psd <- (abs(fft_positive)^2) / N
  
  # Frequency vector for positive frequencies
  frequencies <- seq(0, sampling_freq / 2, length.out = N_half)
  
  # --- Peak Detection within frequency range ---
  freq_idx <- which(frequencies >= freq_range[1] & frequencies <= freq_range[2])
  
  if (length(freq_idx) < 2) {
    stop("Frequency range too narrow or outside signal bandwidth")
  }
  
  # Find maximum peak in the specified frequency range
  psd_range <- psd[freq_idx]
  peak_idx_local <- which.max(psd_range)
  peak_idx_global <- freq_idx[peak_idx_local]
  
  peak_frequency <- frequencies[peak_idx_global]
  peak_power <- psd[peak_idx_global]
  
  # --- Output ---
  result <- list(
    step_freq_hz = peak_frequency,
    step_freq_bpm = peak_frequency * 60,
    peak_frequency = peak_frequency,
    peak_power = peak_power,
    method = "FFT",
    freq_range = freq_range
  )
  
  if (return_spectrum) {
    result$psd <- psd
    result$frequencies <- frequencies
  }
  
  return(result)
}


# =============================================================================
# CONVENIENCE WRAPPER
# =============================================================================

#' Compute Step Frequency (Wrapper)
#'
#' Wrapper function that calls step_frequency_fft() with default parameters.
#'
#' @param acceleration Numeric vector of acceleration signal
#' @param sampling_freq Sampling frequency in Hz (default: 256)
#' @param ... Additional arguments passed to step_frequency_fft()
#'
#' @return List with step frequency results (see step_frequency_fft)
#'
#' @export
step_frequency <- function(acceleration, sampling_freq = 256, ...) {
  step_frequency_fft(acceleration, sampling_freq, ...)
}


# =============================================================================
# VISUALIZATION
# =============================================================================

#' Plot Step Frequency Analysis Results
#'
#' Creates diagnostic plots for step frequency detection showing the time series
#' and the power spectral density with detected peak.
#'
#' @param signal Numeric vector of acceleration signal
#' @param result Result list from step_frequency_fft() with return_spectrum=TRUE
#' @param sampling_freq Sampling frequency in Hz (default: 256)
#' @param xlim_psd X-axis limits for PSD plot in Hz (default: c(0, 5))
#'
#' @return NULL (produces plot as side effect)
#'
#' @examples
#' # Generate test signal
#' t <- seq(0, 10, 1/256)
#' signal <- sin(2 * pi * 1.8 * t) + 0.2 * rnorm(length(t))
#' 
#' # Compute with spectrum
#' result <- step_frequency_fft(signal, return_spectrum = TRUE)
#' 
#' # Plot
#' plot_step_frequency(signal, result)
#'
#' @export
plot_step_frequency <- function(signal, result, sampling_freq = 256, xlim_psd = c(0, 5)) {
  
  # Check if spectrum data is available
  if (is.null(result$psd) || is.null(result$frequencies)) {
    stop("Spectrum data not available. Run step_frequency_fft() with return_spectrum = TRUE")
  }
  
  # Save original par settings and restore on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Set up 2-panel layout with RStudio-safe margins to avoid "margins too large" error
  
  par(mfrow = c(2, 1), 
      mar = c(2.5, 2.5, 1.5, 0.5),
      mgp = c(1.5, 0.4, 0),
      cex = 0.8,
      oma = c(0, 0, 0, 0))
  
  # --- Plot 1: Time series ---
  n_samples <- length(signal)
  time <- seq(0, (n_samples - 1) / sampling_freq, length.out = n_samples)
  
  plot(time, signal, 
       type = "l", 
       col = "steelblue",
       xlab = "Time (s)", 
       ylab = "Acceleration (g)",
       main = "Acceleration Signal")
  grid(col = "gray80")
  
  # Add annotation with detected frequency
  legend("topright", 
         legend = sprintf("%.2f Hz (%.0f steps/min)", 
                          result$step_freq_hz, result$step_freq_bpm),
         text.col = "darkred", 
         bty = "n",
         cex = 0.9)
  
  # --- Plot 2: Power spectral density ---
  # Ensure xlim is within data range
  xlim_actual <- c(
    max(xlim_psd[1], min(result$frequencies)),
    min(xlim_psd[2], max(result$frequencies))
  )
  
  # Find indices within xlim for y-axis scaling
  plot_idx <- which(result$frequencies >= xlim_actual[1] & 
                      result$frequencies <= xlim_actual[2])
  
  if (length(plot_idx) == 0) {
    plot.new()
    text(0.5, 0.5, "No data in specified frequency range", cex = 1.0)
    return(invisible(NULL))
  }
  
  # Scale y-axis based on the search range (where peaks are), not full display range
  search_idx <- which(result$frequencies >= result$freq_range[1] & 
                        result$frequencies <= result$freq_range[2])
  y_max <- max(result$psd[search_idx], na.rm = TRUE)
  
  plot(result$frequencies, result$psd, 
       type = "l", 
       col = "darkblue",
       xlab = "Frequency (Hz)", 
       ylab = "Power",
       main = "Power Spectral Density",
       xlim = xlim_actual,
       ylim = c(0, y_max * 1.1))
  
  # Add grid
  grid(col = "gray80")
  
  # Add shaded region for search range
  if (!is.null(result$freq_range)) {
    rect(result$freq_range[1], 0, result$freq_range[2], y_max * 1.1,
         col = rgb(0, 0.5, 0, 0.1), border = NA)
  }
  
  # Mark the detected peak
  points(result$peak_frequency, result$peak_power, 
         col = "red", pch = 19, cex = 1.5)
  
  # Add vertical line at peak
  abline(v = result$peak_frequency, col = "red", lty = 2, lwd = 1.5)
  
  # Label the peak
  text(result$peak_frequency, result$peak_power * 0.85,
       labels = sprintf("%.2f Hz", result$peak_frequency),
       col = "red", font = 2, pos = 4, cex = 0.9)
  
  invisible(NULL)
}


# =============================================================================
# TEST FUNCTION
# =============================================================================

#' Test Step Frequency Detection
#'
#' Runs a simple test with a synthetic signal to verify the implementation.
#'
#' @param plot_result Logical. Whether to plot the results (default: TRUE)
#'
#' @return List with test results
#'
#' @examples
#' test_step_frequency()
#'
#' @export
test_step_frequency <- function(plot_result = TRUE) {
  
  cat("=== Step Frequency Detection Test ===\n\n")
  
  # Test parameters
  sampling_freq <- 256
  duration <- 10  # seconds
  true_freq <- 1.8  # Hz (108 steps/min)
  
  # Generate test signal: sinusoid + noise
  t <- seq(0, duration, 1/sampling_freq)
  signal <- sin(2 * pi * true_freq * t) + 0.3 * rnorm(length(t))
  
  cat(sprintf("Test signal: %.1f Hz sinusoid + noise\n", true_freq))
  cat(sprintf("Duration: %.1f s, Sampling: %d Hz, Samples: %d\n\n", 
              duration, sampling_freq, length(signal)))
  
  # Compute step frequency
  result <- step_frequency_fft(signal, 
                                sampling_freq = sampling_freq, 
                                return_spectrum = TRUE)
  
  # Report results
  cat("Detection Results:\n")
  cat(sprintf("  Detected frequency: %.3f Hz\n", result$step_freq_hz))
  cat(sprintf("  Expected frequency: %.3f Hz\n", true_freq))
  cat(sprintf("  Error: %.4f Hz (%.2f%%)\n", 
              abs(result$step_freq_hz - true_freq),
              100 * abs(result$step_freq_hz - true_freq) / true_freq))
  cat(sprintf("  Cadence: %.1f steps/min\n", result$step_freq_bpm))
  cat(sprintf("  Search range: %.1f - %.1f Hz\n", result$freq_range[1], result$freq_range[2]))
  
  # Check if detection is accurate (within 5%)
  error_pct <- 100 * abs(result$step_freq_hz - true_freq) / true_freq
  if (error_pct < 5) {
    cat("\n[OK] TEST PASSED: Detection within 5% of true frequency\n")
  } else {
    cat("\n[FAIL] TEST FAILED: Detection error exceeds 5%\n")
  }
  
  # Plot if requested
  if (plot_result) {
    plot_step_frequency(signal, result, sampling_freq)
  }
  
  invisible(result)
}


# =============================================================================
# END OF SCRIPT
# =============================================================================
