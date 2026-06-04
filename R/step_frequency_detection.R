# =============================================================================
# Step Frequency Detection - FFT-Based Method with Parabolic Interpolation
# =============================================================================
#
# RACING Project (ORD 2025) - Version 3.0
# Author: Philippe Terrier, HE-Arc Engineering
# Date: March 2026
#
# Purpose: Estimate step frequency from vertical acceleration using FFT
# power spectral density analysis with sub-bin parabolic interpolation.
#
# Method: Standard FFT-based frequency detection with peak refinement
# via parabolic (quadratic) interpolation of the log-magnitude spectrum.
# This provides sub-bin frequency resolution, making the estimate robust
# to small signal differences that would otherwise shift the peak by
# +/- 1 FFT bin (~0.004 Hz), causing large downstream effects.
#
# Why interpolation matters for RACING:
#   Raw FFT bin resolution: df = fs/N ~ 0.004 Hz for typical gait signals
#   Truncation sensitivity: 1 bin SF difference -> ~50 samples -> different
#     resampling ratio -> cumulative phase drift -> decorrelated signals
#   With interpolation: sub-bin precision (~0.0001 Hz) -> truncation matches
#     regardless of which discrete bin the FFT peak falls into.
#
# Reference MATLAB implementation: SF_det_func.m (ACIER study)
#   - MATLAB returns raw FFT bin frequency (no interpolation)
#   - MATLAB uses peakn-based peak detection (SF_det_func source unavailable;
#     standalone SF_det.m uses peakn order 100, range 0-10 Hz)
#   - R uses which.max within a restricted range [1.4, 3.0] Hz where the
#     step frequency is always the dominant peak (stride < 1.4, harmonics > 3.0)
#   - R adds parabolic interpolation for robustness to floating-point
#     differences in tilt correction that can shift the discrete peak by 1 bin
#
# Scientific references:
#   - Smith, J.O. (2011). Spectral Audio Signal Processing. W3K Publishing.
#     Chapter on "Quadratic Interpolation of Spectral Peaks"
#   - Gasior, M. & Gonzalez, J.L. (2004). Improving FFT frequency
#     measurement resolution by parabolic and Gaussian spectrum interpolation.
#     AIP Conference Proceedings, 732(1), 276-285.
#
# =============================================================================

# =============================================================================
# MAIN FUNCTION: FFT-Based Step Frequency Detection
# =============================================================================

#' Compute Step Frequency Using FFT with Parabolic Peak Interpolation
#'
#' Estimates step frequency from vertical acceleration using FFT power spectral
#' density analysis. After finding the discrete peak bin, parabolic interpolation
#' on the log-magnitude spectrum refines the estimate to sub-bin precision.
#'
#' @param acceleration_vertical Numeric vector of vertical acceleration (g units)
#' @param sampling_freq Sampling frequency in Hz (default: 256)
#' @param freq_range Physiological frequency range c(min, max) in Hz.
#'   Default: c(1.4, 3.0) covers normal walking (84-180 steps/min).
#'   Note: truncate_resample() and the batch pipeline may override this
#'   with config-specified values (default c(1.2, 3.0)).
#' @param detrend Logical. Remove mean from signal (default: TRUE)
#' @param zero_pad Logical. Zero-pad to next power of 2 (default: FALSE).
#'   With interpolation, zero-padding is generally unnecessary.
#' @param return_spectrum Logical. Return full PSD for plotting (default: FALSE)
#'
#' @return List containing:
#'   \item{step_freq_hz}{Step frequency in Hz (interpolated, sub-bin precision)}
#'   \item{step_freq_bpm}{Step frequency in beats per minute}
#'   \item{peak_bin_freq}{Raw bin frequency before interpolation (for diagnostics)}
#'   \item{peak_bin_index}{FFT bin index of discrete peak (1-indexed)}
#'   \item{interp_delta}{Interpolation offset from peak bin (fractional bins)}
#'   \item{peak_power}{Power spectral density at peak}
#'   \item{freq_resolution}{FFT frequency resolution df = fs/N (Hz)}
#'   \item{method}{Detection method identifier}
#'   \item{freq_range}{Search range used}
#'   \item{quality}{List with quality indicators and warnings}
#'   \item{psd}{Power spectral density vector (if return_spectrum=TRUE)}
#'   \item{frequencies}{Frequency vector (if return_spectrum=TRUE)}
#'
#' @details
#' Algorithm steps:
#' 1. Preprocessing: remove mean (detrending)
#' 2. FFT computation (no windowing, matching MATLAB SF_det_func.m)
#' 3. Power spectral density: |X(f)|^2 / N
#' 4. Find discrete peak bin within freq_range
#' 5. Parabolic interpolation on log-magnitude for sub-bin precision:
#'      alpha = log(|X[k-1]|)
#'      beta  = log(|X[k]|)
#'      gamma = log(|X[k+1]|)
#'      delta = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma)
#'      f_interp = (k + delta) * df
#'
#' The log-magnitude interpolation is preferred over magnitude or power
#' interpolation because it gives zero bias for signals analyzed without
#' a window function (rectangular window), which matches the MATLAB code.
#'
#' @examples
#' # Simulate walking signal at 1.8 Hz (108 steps/min)
#' t <- seq(0, 60, 1/256)  # 60 seconds at 256 Hz
#' signal <- sin(2 * pi * 1.8 * t) + 0.2 * rnorm(length(t))
#'
#' # Compute step frequency
#' result <- step_frequency_fft(signal, sampling_freq = 256)
#' cat("Detected frequency:", round(result$step_freq_hz, 4), "Hz\n")
#' cat("Bin frequency:     ", round(result$peak_bin_freq, 4), "Hz\n")
#' cat("Interpolation:     ", round(result$interp_delta, 4), "bins\n")
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

  warnings_list <- character(0)

  if (length(acceleration_vertical) < 4 * sampling_freq) {
    w <- sprintf(
      "Signal length is %.1f seconds. For best accuracy, use signals >= 4 seconds.",
      length(acceleration_vertical) / sampling_freq
    )
    warning(w)
    warnings_list <- c(warnings_list, w)
  }

  # --- Preprocessing ---
  signal <- acceleration_vertical

  # Remove DC component (matching MATLAB: fft(y - mean(y)))
  if (detrend) {
    signal <- signal - mean(signal)
  }

  # --- FFT computation ---
  N <- length(signal)

  # Optional zero padding
  if (zero_pad) {
    N_padded <- 2^ceiling(log2(N))
    signal_padded <- c(signal, rep(0, N_padded - N))
  } else {
    signal_padded <- signal
    N_padded <- N
  }

  # Compute FFT
  fft_result <- fft(signal_padded)

  # Frequency resolution (Hz per bin)
  df <- sampling_freq / N_padded

  # Take only positive frequencies (first half + DC)
  N_half <- floor(N_padded / 2) + 1
  fft_positive <- fft_result[1:N_half]

  # Magnitude spectrum (for interpolation) and power spectrum (for output)
  magnitude <- abs(fft_positive)
  psd <- (magnitude^2) / N_padded

  # Frequency vector: bin k (1-indexed) corresponds to frequency (k-1)*df
  frequencies <- (0:(N_half - 1)) * df

  # --- Peak Detection within frequency range ---
  # Design note: uses which.max (maximum amplitude) rather than MATLAB's
  # peakn + first-peak selection. Within [1.4, 3.0] Hz these are equivalent
  # because stride frequency (~0.9-1.2 Hz) is excluded by the lower bound
  # and step harmonics (~3.6+ Hz) are excluded by the upper bound, making
  # the step frequency always the dominant peak in this range.
  freq_idx <- which(frequencies >= freq_range[1] & frequencies <= freq_range[2])

  if (length(freq_idx) < 3) {
    stop("Frequency range too narrow (need at least 3 bins for interpolation)")
  }

  # Find maximum peak in the specified frequency range
  psd_range <- psd[freq_idx]
  peak_idx_local <- which.max(psd_range)
  peak_idx_global <- freq_idx[peak_idx_local]

  # Raw bin frequency (before interpolation)
  peak_bin_freq <- frequencies[peak_idx_global]

  # --- Parabolic interpolation on log-magnitude ---
  # Refine peak position to sub-bin precision using the three-point
  # parabolic interpolation method on log(|X|).
  #
  # For a rectangular window (no windowing), log-magnitude interpolation
  # gives the best accuracy (Gasior & Gonzalez, 2004).
  #
  # delta = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma)
  # where alpha = log|X[k-1]|, beta = log|X[k]|, gamma = log|X[k+1]|
  #
  # The interpolated frequency is: f = (k-1 + delta) * df (0-indexed)
  # In R 1-indexed: f = (peak_idx_global - 1 + delta) * df

  interp_delta <- 0  # default: no correction
  interp_valid <- FALSE

  # Need neighbors on both sides (not at edges of spectrum)
  if (peak_idx_global > 1 && peak_idx_global < N_half) {

    mag_left  <- magnitude[peak_idx_global - 1]
    mag_peak  <- magnitude[peak_idx_global]
    mag_right <- magnitude[peak_idx_global + 1]

    # Guard against zero or near-zero magnitudes (log would fail)
    if (mag_left > 0 && mag_peak > 0 && mag_right > 0) {

      alpha <- log(mag_left)
      beta  <- log(mag_peak)
      gamma <- log(mag_right)

      denom <- alpha - 2 * beta + gamma

      # Denominator must be negative (concave peak) and non-zero
      if (denom < -1e-12) {
        interp_delta <- 0.5 * (alpha - gamma) / denom

        # Sanity check: delta should be in (-0.5, +0.5)
        # Values outside this range indicate a non-parabolic peak shape
        if (abs(interp_delta) <= 0.5) {
          interp_valid <- TRUE
        } else {
          # Clamp to valid range and flag warning
          interp_delta <- max(-0.5, min(0.5, interp_delta))
          w <- sprintf(
            "Interpolation delta %.3f outside normal range; clamped to %.1f",
            interp_delta, sign(interp_delta) * 0.5)
          warnings_list <- c(warnings_list, w)
          interp_valid <- TRUE  # still use clamped value
        }
      } else {
        # Flat or convex peak neighborhood: skip interpolation
        w <- "Interpolation skipped: non-concave peak (flat spectral plateau?)"
        warnings_list <- c(warnings_list, w)
      }
    } else {
      w <- "Interpolation skipped: zero magnitude at peak or neighbor"
      warnings_list <- c(warnings_list, w)
    }
  } else {
    w <- "Interpolation skipped: peak at edge of frequency range"
    warnings_list <- c(warnings_list, w)
  }

  # Compute interpolated frequency
  # peak_idx_global is 1-indexed: frequency = (peak_idx_global - 1) * df
  # With interpolation: frequency = (peak_idx_global - 1 + delta) * df
  peak_frequency <- (peak_idx_global - 1 + interp_delta) * df

  # Peak power (at discrete bin)
  peak_power <- psd[peak_idx_global]

  # --- Quality metrics ---
  quality <- list(
    interp_valid = interp_valid,
    interp_delta = interp_delta,
    peak_bin_freq = peak_bin_freq,
    freq_resolution = df,
    warnings = warnings_list
  )

  # --- Output ---
  result <- list(
    step_freq_hz   = peak_frequency,
    step_freq_bpm  = peak_frequency * 60,
    peak_bin_freq  = peak_bin_freq,
    peak_bin_index = peak_idx_global,
    interp_delta   = interp_delta,
    peak_power     = peak_power,
    freq_resolution = df,
    method         = "FFT_parabolic",
    freq_range     = freq_range
  )

  # Attach quality info
  result$quality <- quality

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
#' Creates diagnostic plots showing the time series and the power spectral
#' density with detected peak and interpolation.
#'
#' @param signal Numeric vector of acceleration signal
#' @param result Result list from step_frequency_fft() with return_spectrum=TRUE
#' @param sampling_freq Sampling frequency in Hz (default: 256)
#' @param xlim_psd X-axis limits for PSD plot in Hz (default: c(0, 5))
#'
#' @return NULL (produces plot as side effect)
#'
#' @export
plot_step_frequency <- function(signal, result, sampling_freq = 256,
                                xlim_psd = c(0, 5)) {

  # Check if spectrum data is available
  if (is.null(result$psd) || is.null(result$frequencies)) {
    stop("Spectrum data not available. Run step_frequency_fft() with return_spectrum = TRUE")
  }

  # Save/restore graphics parameters (prevents margin errors)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

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

  legend("topright",
         legend = sprintf("%.4f Hz (%.1f steps/min)",
                          result$step_freq_hz, result$step_freq_bpm),
         text.col = "darkred",
         bty = "n",
         cex = 0.9)

  # --- Plot 2: Power spectral density ---
  xlim_actual <- c(
    max(xlim_psd[1], min(result$frequencies)),
    min(xlim_psd[2], max(result$frequencies))
  )

  plot_idx <- which(result$frequencies >= xlim_actual[1] &
                      result$frequencies <= xlim_actual[2])

  if (length(plot_idx) == 0) {
    plot.new()
    text(0.5, 0.5, "No data in specified frequency range", cex = 1.0)
    return(invisible(NULL))
  }

  # Scale y-axis based on search range
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

  grid(col = "gray80")

  # Shaded search range
  if (!is.null(result$freq_range)) {
    rect(result$freq_range[1], 0, result$freq_range[2], y_max * 1.1,
         col = rgb(0, 0.5, 0, 0.1), border = NA)
  }

  # Mark the discrete peak bin (gray circle)
  bin_idx <- result$peak_bin_index
  points(result$peak_bin_freq, result$psd[bin_idx],
         col = "gray50", pch = 1, cex = 1.8, lwd = 1.5)

  # Mark the interpolated peak (red filled circle)
  abline(v = result$step_freq_hz, col = "red", lty = 2, lwd = 1.5)
  points(result$step_freq_hz, result$peak_power,
         col = "red", pch = 19, cex = 1.3)

  # Label
  label_text <- sprintf("%.4f Hz (bin: %.4f, delta: %+.3f)",
                        result$step_freq_hz,
                        result$peak_bin_freq,
                        result$interp_delta)
  text(result$step_freq_hz, result$peak_power * 0.80,
       labels = label_text,
       col = "red", font = 2, pos = 4, cex = 0.85)

  invisible(NULL)
}


# =============================================================================
# TEST FUNCTION
# =============================================================================

#' Test Step Frequency Detection with Interpolation
#'
#' Tests the implementation with a synthetic signal at a known frequency
#' that falls between two FFT bins, verifying interpolation accuracy.
#'
#' @param plot_result Logical. Whether to plot the results (default: TRUE)
#'
#' @return List with test results
#'
#' @export
test_step_frequency <- function(plot_result = TRUE) {

  cat("=== Step Frequency Detection Test (with interpolation) ===\n\n")

  # Test parameters
  sampling_freq <- 256
  duration <- 60  # seconds (long signal = fine bin resolution)
  true_freq <- 1.8765  # Hz - deliberately between bins

  # Generate test signal: sinusoid + harmonics + noise (gait-like)
  t <- seq(0, duration, 1/sampling_freq)
  signal <- sin(2 * pi * true_freq * t) +
            0.4 * sin(2 * pi * 2 * true_freq * t) +  # stride harmonic
            0.2 * rnorm(length(t))

  N <- length(signal)
  df <- sampling_freq / N
  cat(sprintf("Test signal: %.4f Hz + stride harmonic + noise\n", true_freq))
  cat(sprintf("Duration: %.1f s, N = %d, df = %.6f Hz\n", duration, N, df))
  cat(sprintf("Nearest bins: %.6f Hz and %.6f Hz\n\n",
              floor(true_freq / df) * df, ceiling(true_freq / df) * df))

  # Compute step frequency
  result <- step_frequency_fft(signal,
                                sampling_freq = sampling_freq,
                                return_spectrum = TRUE)

  # Report results
  cat("Detection Results:\n")
  cat(sprintf("  Bin frequency:    %.6f Hz (bin %d)\n",
              result$peak_bin_freq, result$peak_bin_index))
  cat(sprintf("  Interp delta:     %+.4f bins\n", result$interp_delta))
  cat(sprintf("  Interp frequency: %.6f Hz\n", result$step_freq_hz))
  cat(sprintf("  True frequency:   %.6f Hz\n", true_freq))
  cat(sprintf("  Bin error:        %.4f Hz (%.3f bins)\n",
              abs(result$peak_bin_freq - true_freq),
              abs(result$peak_bin_freq - true_freq) / df))
  cat(sprintf("  Interp error:     %.4f Hz (%.3f bins)\n",
              abs(result$step_freq_hz - true_freq),
              abs(result$step_freq_hz - true_freq) / df))
  cat(sprintf("  Improvement:      %.1fx\n",
              abs(result$peak_bin_freq - true_freq) /
              max(abs(result$step_freq_hz - true_freq), 1e-10)))

  # Check truncation equivalence
  # Note: the actual pipeline rounds SF to 0.1 Hz before computing
  # the truncation index (see truncate_resample.R). This absorbs small
  # bin differences and ensures identical MATLAB/R truncation lengths.
  trunc_true <- floor(256 / true_freq * 250)
  trunc_bin  <- floor(256 / result$peak_bin_freq * 250)
  trunc_int  <- floor(256 / result$step_freq_hz * 250)
  trunc_rnd  <- floor(256 / round(result$step_freq_hz, 1) * 250)
  cat(sprintf("\n  Truncation (true freq):   %d samples\n", trunc_true))
  cat(sprintf("  Truncation (bin freq):    %d samples (diff: %d)\n",
              trunc_bin, trunc_bin - trunc_true))
  cat(sprintf("  Truncation (interp freq): %d samples (diff: %d)\n",
              trunc_int, trunc_int - trunc_true))
  cat(sprintf("  Truncation (rounded 0.1): %d samples (pipeline uses this)\n",
              trunc_rnd))

  # Test: verify that two signals with +1 bin shift converge
  cat("\n--- Bin boundary robustness test ---\n")
  # Shift signal slightly to make peak fall on the other side of boundary
  signal2 <- sin(2 * pi * true_freq * t + 0.01) +  # tiny phase shift
             0.4 * sin(2 * pi * 2 * true_freq * t + 0.01) +
             0.2 * rnorm(length(t))

  result2 <- step_frequency_fft(signal2, sampling_freq = sampling_freq)

  cat(sprintf("  Signal 1: bin %d (%.6f Hz) -> interp %.6f Hz\n",
              result$peak_bin_index, result$peak_bin_freq, result$step_freq_hz))
  cat(sprintf("  Signal 2: bin %d (%.6f Hz) -> interp %.6f Hz\n",
              result2$peak_bin_index, result2$peak_bin_freq, result2$step_freq_hz))
  cat(sprintf("  Bin difference: %d\n",
              result2$peak_bin_index - result$peak_bin_index))
  cat(sprintf("  Interp difference: %.6f Hz (%.3f bins)\n",
              abs(result2$step_freq_hz - result$step_freq_hz),
              abs(result2$step_freq_hz - result$step_freq_hz) / df))

  trunc1 <- floor(256 / result$step_freq_hz * 250)
  trunc2 <- floor(256 / result2$step_freq_hz * 250)
  trunc1_r <- floor(256 / round(result$step_freq_hz, 1) * 250)
  trunc2_r <- floor(256 / round(result2$step_freq_hz, 1) * 250)
  cat(sprintf("  Truncation match:  %s (diff: %d samples)\n",
              if (trunc1 == trunc2) "[OK]" else "[DIFF]",
              abs(trunc1 - trunc2)))
  cat(sprintf("  Rounded match:     %s (diff: %d, pipeline uses this)\n",
              if (trunc1_r == trunc2_r) "[OK]" else "[DIFF]",
              abs(trunc1_r - trunc2_r)))

  # Plot if requested
  if (plot_result) {
    plot_step_frequency(signal, result, sampling_freq)
  }

  invisible(result)
}


# =============================================================================
# END OF SCRIPT
# =============================================================================
