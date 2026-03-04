#' ============================================================================
#' Gait Signal Truncation and Resampling Pipeline
#' ============================================================================
#' 
#' Modular preprocessing function for accelerometer gait signals. Applies tilt
#' correction, truncates signals to a specified number of steps, and resamples 
#' to a standardized length for nonlinear dynamics analysis.
#' 
#' @description
#' This function produces two output versions of each signal:
#' 1. TRUNCATED: Original sampling rate, containing exactly N steps
#'    - Used for: RMS, step frequency, ACF regularity (temporal metrics)
#'    
#' 2. RESAMPLED: Standardized length with uniform samples per step
#'    - Used for: AMI, LDS, ACI (phase space reconstruction metrics)
#' 
#' The dual-output design is necessary because linear metrics depend on actual
#' temporal relationships, while nonlinear metrics require standardized time
#' scales for valid inter-subject comparisons.
#' 
#' ## SCIENTIFIC REFERENCES
#'
#' Tilt correction algorithm (Moe-Nilssen method):
#'   Moe-Nilssen, R. (1998). A new method for evaluating motor control in gait 
#'   under real-life environmental conditions. Part 1: The instrument.
#'   Clinical Biomechanics, 13(4-5), 320-327.
#'   https://doi.org/10.1016/S0268-0033(98)00089-8
#'
#' ACIER study methodology and validation:
#'   Piergiovanni, S., & Terrier, P. (2024).
#'   Effects of metronome walking on long-term attractor divergence and 
#'   correlation structure of gait: A validation study in older people. 
#'   Scientific Reports, 14, 15784.
#'   https://doi.org/10.1038/s41598-024-65662-5
#'
#'   Piergiovanni, S., & Terrier, P. (2024).
#'   Validity of linear and nonlinear measures of gait variability to 
#'   characterize aging gait with a single lower back accelerometer.
#'   Sensors, 24(23), 7427.
#'   https://doi.org/10.3390/s24237427
#'
#' Signal resampling methodology:
#'   Smith, J. O., & Gossett, P. (1984).
#'   A flexible sampling-rate conversion method. In Proceedings of ICASSP '84.
#'   IEEE International Conference on Acoustics, Speech, and Signal Processing.
#'
#' @author Philippe Terrier / RACING Project (ORD 2025), HE-Arc 
#' @date February 2026
#' ============================================================================

# =============================================================================
# MAIN FUNCTION: Truncate and Resample Gait Signal
# =============================================================================

#' Truncate and Resample Gait Acceleration Signals
#' 
#' @description
#' Preprocesses triaxial accelerometer data by:
#' 1. Applying tilt correction (Moe-Nilssen 1998) to remove gravity effects
#' 2. Detecting step frequency from vertical acceleration
#' 3. Truncating to a specified number of steps
#' 4. Computing the acceleration norm (vector magnitude)
#' 5. Resampling to a standardized length
#' 
#' If the input signal is too short to contain the requested number of steps,
#' the function uses all available samples and issues a warning.
#' 
#' @param acc_x Numeric vector. Anteroposterior (AP) acceleration in g units.
#'   Raw signal (tilt correction will be applied if apply_tilt_correction=TRUE).
#' @param acc_y Numeric vector. Vertical (V) acceleration in g units.
#'   Raw signal. Used for step frequency detection after tilt correction.
#' @param acc_z Numeric vector. Mediolateral (ML) acceleration in g units.
#'   Raw signal (tilt correction will be applied if apply_tilt_correction=TRUE).
#' @param sampling_freq Numeric. Sampling frequency in Hz. Default: 256.
#' @param target_steps Integer. Number of steps to extract. Default: 250.
#'   ACIER Common values: 250 (NC/AC conditions), 500 (outdoor condition).
#' @param target_samples Integer or NULL. Target length after resampling.
#'   If NULL (default), computed as target_steps * samples_per_step.
#' @param samples_per_step Integer. Number of samples per step after resampling.
#'   Default: 75 (ACIER standard). Used only if target_samples is NULL.
#' @param apply_tilt_correction Logical. If TRUE (default), apply Moe-Nilssen
#'   (1998) tilt correction before processing. Set FALSE if signals are already
#'   tilt-corrected.
#' @param auto_orient Logical. If TRUE (default) and tilt correction is applied,
#'   automatically detect and correct for inverted sensor orientation.
#' @param step_freq_hz Numeric or NULL. Pre-computed step frequency in Hz.
#'   When apply_tilt_correction=TRUE, this value is IGNORED and SF is always
#'   re-detected internally on the tilt-corrected vertical signal, matching
#'   the MATLAB/ACIER pipeline order (calibphys -> SF_det_func). This prevents
#'   a 1-FFT-bin shift (~0.003 Hz) that cascades into resampling phase drift.
#'   When apply_tilt_correction=FALSE (pre-corrected data), this value is
#'   used directly, skipping internal detection.
#' @param freq_range Numeric vector of length 2. Physiological frequency range
#'   c(min, max) in Hz for step frequency detection. Default: c(1.2, 3.0) 
#'   covers typical adult gait (72-180 steps/min). Only ignored when
#'   step_freq_hz is provided AND apply_tilt_correction=FALSE.
#' @param compute_norm Logical. If TRUE (default), compute and include 
#'   the vector norm (magnitude) in outputs.
#' @param norm_subtract_gravity Logical. If TRUE (default), subtract 1g from
#'   the norm (sqrt(x^2 + y^2 + z^2) - 1). Set FALSE if you want the raw magnitude.
#' @param resampling_method Character. Method for resampling: "auto" (default),
#'   "direct", or "two_stage". See resample_signal() for details.
#' @param verbose Logical. If TRUE, print progress messages. Default: FALSE.
#' 
#' @return A list with the following components:
#' \describe{
#'   \item{truncated}{List containing truncated signals at original sampling rate:
#'     \itemize{
#'       \item \code{x}: AP acceleration (tilt-corrected if applied)
#'       \item \code{y}: Vertical acceleration (tilt-corrected if applied)
#'       \item \code{z}: ML acceleration (tilt-corrected if applied)
#'       \item \code{norm}: Vector magnitude (if compute_norm = TRUE)
#'     }
#'   }
#'   \item{resampled}{List containing resampled signals at standardized length:
#'     \itemize{
#'       \item \code{x}: AP acceleration
#'       \item \code{y}: Vertical acceleration
#'       \item \code{z}: ML acceleration
#'       \item \code{norm}: Vector magnitude (if compute_norm = TRUE)
#'     }
#'   }
#'   \item{metadata}{List of processing parameters:
#'     \itemize{
#'       \item \code{step_freq_hz}: Detected or provided step frequency
#'       \item \code{step_freq_bpm}: Step frequency in beats per minute
#'       \item \code{target_steps}: Requested number of steps
#'       \item \code{actual_steps}: Actual number of steps in truncated signal
#'       \item \code{truncated_length}: Length of truncated signals
#'       \item \code{resampled_length}: Length of resampled signals
#'       \item \code{samples_per_step}: Samples per step after resampling
#'       \item \code{sampling_freq}: Original sampling frequency
#'       \item \code{resampling_ratio}: Ratio target_samples / truncated_length
#'       \item \code{resampling_method}: Method used ("direct" or "two_stage")
#'       \item \code{signal_truncated}: Was the signal actually truncated?
#'       \item \code{signal_too_short}: Was the signal shorter than requested?
#'       \item \code{tilt_correction_applied}: Was tilt correction applied?
#'     }
#'   }
#'   \item{tilt_correction}{List of tilt correction diagnostics (if applied):
#'     \itemize{
#'       \item \code{theta_ap_deg}: AP tilt angle in degrees
#'       \item \code{theta_ml_deg}: ML tilt angle in degrees
#'       \item \code{orientation_detected}: "upward" or "downward"
#'       \item \code{orientation_corrected}: Was orientation flipped?
#'     }
#'   }
#'   \item{quality}{List of quality indicators:
#'     \itemize{
#'       \item \code{tilt_correction_valid}: Did tilt correction pass validation?
#'       \item \code{warnings}: Character vector of any warnings generated
#'     }
#'   }
#' }
#' 
#' 
#' @details
#' ## Processing Pipeline (Exact ACIER/MATLAB Order)
#' 
#' ```
#' Raw Signals (x, y, z)
#'         |
#'         v
#' [0] Compute Norm from RAW signals -----> norm_raw
#'         |
#'         v
#' [1] Tilt Correction (Moe-Nilssen) ----> (x', y', z')
#'         |
#'         v
#' [2] Step Frequency Detection (FFT) ----> SF [Hz]
#'     ** MUST use tilt-corrected y' **
#'     ** (even if step_freq_hz provided externally) **
#'         |
#'         v
#' [3] Calculate Truncation Index --------> N_trunc
#'         |
#'         v
#' [4] Truncate All Signals --------------> truncated output
#'         |
#'         v
#' [5] Resample to Standard Length -------> resampled output
#' ```
#' 
#' NOTE: The order [1] before [2] is critical. MATLAB applies tilt correction
#' (calibphys) before SF detection (SF_det_func). Detecting SF on raw
#' (uncorrected) data shifts the FFT peak by 1 bin (~0.003 Hz) in ~50% of
#' subjects, causing ~50 samples truncation difference, different resampling
#' ratios, cumulative phase drift (r=0.99 -> r=0.02), and AMI/LDS/ACI
#' mismatches.
#' 
#' ## Signal Length Calculations
#' 
#' 1. **Tilt Correction** (optional): Applies Moe-Nilssen (1998) algorithm to
#'    remove the effects of sensor tilt from gravity. Automatically detects
#'    inverted sensor orientation.
#' 
#' 2. **Step Frequency Detection**: Uses FFT on vertical acceleration to find
#'    the dominant walking frequency. When tilt correction is applied, SF is
#'    always detected internally on the corrected signal (even if step_freq_hz
#'    is provided). Can only be skipped when apply_tilt_correction=FALSE.
#' 
#' 3. **Truncation Index Calculation**: 
#'    \code{truncation_index = floor(sampling_freq / step_freq_hz * target_steps)}
#' 
#' 4. **Short Signal Handling**: If the signal is shorter than truncation_index,
#'    the function uses all available samples. The actual number of steps
#'    is computed and reported in metadata.
#' 
#' 5. **Norm Computation**: Vector magnitude is computed as:
#'    - With gravity subtraction: sqrt(x^2 + y^2 + z^2) - 1
#'    - Without gravity subtraction: sqrt(x^2 + y^2 + z^2)
#' 
#' 6. **Resampling**: Uses bandlimited interpolation (signal::resample) to
#'    achieve the target length while preserving signal dynamics. Automatically
#'    handles problematic ratios close to 1.0 using two-stage resampling.
#' 
#' @examples
#' \dontrun{
#' # Load accelerometer data (assumes CSV with columns: acc_x, acc_y, acc_z)
#' data <- read.csv("subject_001_NC.csv")
#' 
#' # ==================================================================
#' # Example 1: Standard ACIER NC condition (250 steps)
#' # ==================================================================
#' result_nc <- truncate_resample(
#'   acc_x = data$acc_x,
#'   acc_y = data$acc_y,
#'   acc_z = data$acc_z,
#'   sampling_freq = 256,
#'   target_steps = 250,
#'   samples_per_step = 75,        # ACIER standard
#'   apply_tilt_correction = TRUE,
#'   verbose = TRUE
#' )
#' 
#' # Access outputs
#' print(result_nc)                           # Summary
#' truncated_ml <- result_nc$truncated$z      # For RMS, ACF
#' resampled_norm <- result_nc$resampled$norm # For AMI, LDS, ACI
#' sf_hz <- result_nc$metadata$step_freq_hz   # Detected frequency
#' 
#' # ==================================================================
#' # Example 2: ACIER Outdoor condition (500 steps)
#' # ==================================================================
#' result_out <- truncate_resample(
#'   acc_x = data$acc_x,
#'   acc_y = data$acc_y,
#'   acc_z = data$acc_z,
#'   target_steps = 500,
#'   target_samples = 37500,  # 500  x  75
#'   verbose = TRUE
#' )
#' 
#' # ==================================================================
#' # Example 3: Pre-corrected data with known step frequency
#' # ==================================================================
#' result_fast <- truncate_resample(
#'   acc_x = data$acc_x_corrected,
#'   acc_y = data$acc_y_corrected,
#'   acc_z = data$acc_z_corrected,
#'   apply_tilt_correction = FALSE,  # Already corrected
#'   step_freq_hz = 1.85,             # Skip SF detection
#'   target_steps = 250,
#'   verbose = FALSE
#' )
#' 
#' # ==================================================================
#' # Example 4: Quality checking
#' # ==================================================================
#' if (length(result_nc$quality$warnings) > 0) {
#'   cat("Warnings detected:\n")
#'   print(result_nc$quality$warnings)
#' }
#' 
#' if (!result_nc$quality$tilt_correction_valid) {
#'   warning("Tilt correction quality check failed - review data")
#' }
#' 
#' # Check actual vs target steps
#' if (result_nc$metadata$signal_too_short) {
#'   cat(sprintf("Only %.1f steps available (target: %d)\n",
#'               result_nc$metadata$actual_steps,
#'               result_nc$metadata$target_steps))
#' }
#' 
#' # ==================================================================
#' # Example 5: Batch processing with error handling
#' # ==================================================================
#' process_subject <- function(subject_id, condition) {
#'   tryCatch({
#'     data <- read.csv(sprintf("sub_%03d_%s.csv", subject_id, condition))
#'     
#'     result <- truncate_resample(
#'       acc_x = data$acc_x,
#'       acc_y = data$acc_y,
#'       acc_z = data$acc_z,
#'       target_steps = 250,
#'       verbose = FALSE
#'     )
#'     
#'     # Extract key metrics for downstream analysis
#'     list(
#'       subject_id = subject_id,
#'       condition = condition,
#'       step_freq_hz = result$metadata$step_freq_hz,
#'       resampled_norm = result$resampled$norm,
#'       quality_ok = result$quality$tilt_correction_valid &&
#'                    !result$metadata$signal_too_short
#'     )
#'   }, error = function(e) {
#'     warning(sprintf("Failed for subject %d, condition %s: %s",
#'                     subject_id, condition, e$message))
#'     NULL
#'   })
#' }
#' 
#' # Process multiple subjects
#' subjects <- 1:60
#' conditions <- c("NC", "AC")
#' results <- lapply(subjects, function(id) {
#'   lapply(conditions, function(cond) {
#'     process_subject(id, cond)
#'   })
#' })
#' }
#' 
#' @seealso 
#' \code{\link{tilt_correction}} for the Moe-Nilssen tilt correction algorithm
#' \code{\link{resample_signal}} for the resampling implementation
#' \code{\link{step_frequency_fft}} for FFT-based step frequency detection
#' 
#' @export
truncate_resample <- function(acc_x, 
                               acc_y, 
                               acc_z,
                               sampling_freq = 256,
                               target_steps = 250,
                               target_samples = NULL,
                               samples_per_step = 75,
                               apply_tilt_correction = TRUE,
                               auto_orient = TRUE,
                               step_freq_hz = NULL,
                               freq_range = c(1.2, 3.0), 
                               compute_norm = TRUE,
                               norm_subtract_gravity = TRUE,
                               resampling_method = c("auto", "direct", "two_stage"),
                               verbose = TRUE,
                               log_file = NULL) {
  
  # ===========================================================================
  # Input Validation
  # ===========================================================================
  
  resampling_method <- match.arg(resampling_method)
  warnings_list <- character(0)
  
  # Check that all inputs are numeric vectors of the same length
  if (!is.numeric(acc_x) || !is.numeric(acc_y) || !is.numeric(acc_z)) {
    stop("All acceleration inputs must be numeric vectors")
  }
  
  n <- length(acc_x)
  if (length(acc_y) != n || length(acc_z) != n) {
    stop("All acceleration vectors must have the same length")
  }
  
  if (n < 256) {
    stop("Signal too short: minimum 256 samples (1 second at 256 Hz) required")
  }
  
  # Validate parameters
  if (!is.numeric(sampling_freq) || sampling_freq <= 0) {
    stop("sampling_freq must be a positive number")
  }
  
  if (!is.numeric(target_steps) || target_steps < 1) {
    stop("target_steps must be a positive integer")
  }
  target_steps <- as.integer(target_steps)
  
  if (!is.null(target_samples)) {
    if (!is.numeric(target_samples) || target_samples < 1) {
      stop("target_samples must be a positive integer or NULL")
    }
    target_samples <- as.integer(target_samples)
  }
  
  if (!is.numeric(samples_per_step) || samples_per_step < 1) {
    stop("samples_per_step must be a positive integer")
  }
  samples_per_step <- as.integer(samples_per_step)
  
  # Check for NA values
  if (any(is.na(acc_x)) || any(is.na(acc_y)) || any(is.na(acc_z))) {
    warning("NA values detected in input signals. These may cause issues.")
    warnings_list <- c(warnings_list, "NA values in input signals")
  }
  
  if (verbose) {
    message(sprintf("Input signal: %d samples (%.1f seconds at %d Hz)",
                    n, n / sampling_freq, sampling_freq))
  }
  
  # ===========================================================================
  # Step 0: Compute Norm BEFORE Tilt Correction
  # ===========================================================================
  # IMPORTANT: Following MATLAB/ACIER implementation, the vector norm must be
  # computed on RAW signals BEFORE tilt correction.
  #
  # MATHEMATICAL NOTE: Tilt correction applies a rotation matrix R to align
  # sensor axes with anatomical directions. Since rotations preserve vector
  # magnitude (||R*v|| = ||v||), the norm is mathematically invariant to tilt
  # correction. However, we compute it BEFORE correction to exactly replicate
  # the MATLAB processing order for validation purposes.
  #
  # ACIER VALIDATION REQUIREMENT: For R-MATLAB comparison studies, this order
  # must be preserved. Computing norm after tilt correction would give identical
  # results but break byte-level reproducibility with MATLAB outputs.
  
  norm_full_raw <- NULL
  if (compute_norm) {
    if (verbose) message("Computing vector norm (before tilt correction)...")
    norm_full_raw <- sqrt(acc_x^2 + acc_y^2 + acc_z^2)
    if (norm_subtract_gravity) {
      norm_full_raw <- norm_full_raw - 1
    }
  }
  
  # ===========================================================================
  # Step 1: Tilt Correction (Moe-Nilssen 1998)
  # ===========================================================================
  
  tilt_diagnostics <- NULL
  tilt_correction_valid <- NA
  
  if (apply_tilt_correction) {
    if (verbose) message("Applying tilt correction (Moe-Nilssen 1998)...")
    
    # Check if tilt_correction function is available
    if (!exists("tilt_correction", mode = "function")) {
      # Try to source it from the current directory
      if (file.exists("tilt_correction.R")) {
        source("tilt_correction.R")
      } else {
        stop("tilt_correction() function not found. Source tilt_correction.R first.")
      }
    }
    
    # Apply tilt correction with diagnostics
    tilt_result <- tilt_correction(
      acc_ap = acc_x,
      acc_v = acc_y,
      acc_ml = acc_z,
      auto_orient = auto_orient,
      return_diagnostics = TRUE,
      verbose = verbose
    )
    
    # Update acceleration signals with corrected values
    acc_x <- tilt_result$ap
    acc_y <- tilt_result$v
    acc_z <- tilt_result$ml
    
    # Store diagnostics
    tilt_diagnostics <- list(
      theta_ap_deg = tilt_result$theta_ap_deg,
      theta_ml_deg = tilt_result$theta_ml_deg,
      orientation_detected = tilt_result$orientation_detected,
      orientation_corrected = tilt_result$orientation_corrected
    )
    
    tilt_correction_valid <- tilt_result$quality_check
    
    if (verbose) {
      message(sprintf("  Tilt angles: AP = %.1f deg, ML = %.1f deg",
                      tilt_result$theta_ap_deg, tilt_result$theta_ml_deg))
      if (tilt_result$orientation_corrected) {
        message("  Note: Inverted sensor orientation was auto-corrected")
      }
    }
    
    # Quality warnings from tilt correction
    if (!tilt_correction_valid) {
      warn_msg <- sprintf("Tilt correction quality check failed (output means not near zero)")
      warnings_list <- c(warnings_list, warn_msg)
      if (verbose) message(sprintf("  WARNING: %s", warn_msg))
    }
    
    # Extreme tilt warning
    if (abs(tilt_result$theta_ap_deg) > 30 || abs(tilt_result$theta_ml_deg) > 30) {
      warn_msg <- sprintf("Extreme tilt detected: AP=%.1f deg, ML=%.1f deg",
                          tilt_result$theta_ap_deg, tilt_result$theta_ml_deg)
      warnings_list <- c(warnings_list, warn_msg)
      if (verbose) message(sprintf("  WARNING: %s", warn_msg))
    }
    
  } else {
    if (verbose) message("Skipping tilt correction (already corrected or not needed)")
  }
  
  # ===========================================================================
  # Step 2: Step Frequency Detection
  # ===========================================================================
  # CRITICAL: SF must be detected on the TILT-CORRECTED vertical signal to
  # match the MATLAB/ACIER pipeline (calibphys -> SF_det_func). Tilt
  # correction redistributes energy between axes, shifting the FFT peak by
  # up to 1 frequency bin (~0.003 Hz). This cascades into truncation length
  # differences (~50 samples), different resampling ratios, and progressive
  # phase drift that decorrelates the resampled signals (r: 0.99 -> 0.02).
  #
  # Therefore, when apply_tilt_correction=TRUE, SF is ALWAYS detected
  # internally on the corrected acc_y, even if step_freq_hz was provided
  # externally (which may have been computed on uncorrected data).
  #
  # When apply_tilt_correction=FALSE (pre-corrected data), an externally
  # provided step_freq_hz is used as-is.
  # ===========================================================================
  
  # Determine whether SF must be (re-)detected internally
  detect_sf_internally <- is.null(step_freq_hz) ||
                          (apply_tilt_correction && !is.null(step_freq_hz))
  
  if (detect_sf_internally) {
    
    if (!is.null(step_freq_hz) && apply_tilt_correction && verbose) {
      message(sprintf(
        "  Note: Overriding provided SF (%.4f Hz) - re-detecting on tilt-corrected signal",
        step_freq_hz))
    }
    
    if (verbose) message("Detecting step frequency from tilt-corrected vertical acceleration...")
    
    # Check if step_frequency_fft function is available
    if (!exists("step_frequency_fft", mode = "function")) {
      # Try to source it from the current directory
      if (file.exists("step_frequency_detection.R")) {
        source("step_frequency_detection.R")
      } else {
        stop("step_frequency_fft() function not found. Source step_frequency_detection.R first.")
      }
    }
    
    sf_result <- step_frequency_fft(acc_y, 
                                    sampling_freq = sampling_freq,
                                    freq_range = freq_range)
    
    step_freq_hz <- sf_result$step_freq_hz
    
    if (verbose) {
      message(sprintf("  Detected: %.3f Hz (%.1f steps/min)",
                      step_freq_hz, step_freq_hz * 60))
    }
  } else {
    if (verbose) {
      message(sprintf("  Using provided SF: %.4f Hz (%.1f steps/min) [no tilt correction]",
                      step_freq_hz, step_freq_hz * 60))
    }
  }
  
  step_freq_bpm <- step_freq_hz * 60
  
  # Round SF to 0.1 Hz for truncation (RACING: ensures MATLAB/R convergence)
  # Precise SF is preserved in metadata for reporting
  step_freq_hz_precise <- step_freq_hz
  step_freq_hz <- round(step_freq_hz, 1)
  
  # Log step frequency to debug file (if active)
  resample_log(log_file, "SF detection: precise=%.4f Hz, rounded=%.1f Hz (%.1f bpm)",
               step_freq_hz_precise, step_freq_hz, step_freq_bpm)
  
  # ===========================================================================
  # Step 3: Calculate Truncation Index
  # ===========================================================================
  
  # Number of samples per step at original sampling rate
  samples_per_step_original <- sampling_freq / step_freq_hz
  
  # Total samples needed for target_steps
  truncation_index <- floor(samples_per_step_original * target_steps)
  
  # Log truncation to debug file (if active)
  resample_log(log_file, "Truncation: %.2f samples/step x %d steps = %d samples (signal has %d)",
               samples_per_step_original, target_steps, truncation_index, n)
  
  if (verbose) {
    message(sprintf("Truncation: %.1f samples/step  x  %d steps = %d samples needed",
                    samples_per_step_original, target_steps, truncation_index))
  }
  
  # ===========================================================================
  # Step 4: Handle Short Signals
  # ===========================================================================
  
  signal_too_short <- FALSE
  signal_truncated <- TRUE
  actual_steps <- target_steps
  
  if (truncation_index > n) {
    # Signal is too short - use all available samples
    signal_too_short <- TRUE
    signal_truncated <- FALSE
    truncation_index <- n
    actual_steps <- n / samples_per_step_original
    
    warning_msg <- sprintf(
      "Signal too short for %d steps (need %d samples, have %d). Using all %d samples (%.1f steps).",
      target_steps, 
      floor(samples_per_step_original * target_steps),
      n,
      n,
      actual_steps
    )
    warning(warning_msg)
    warnings_list <- c(warnings_list, warning_msg)
    
    if (verbose) {
      message(sprintf("  WARNING: %s", warning_msg))
    }
  } else if (truncation_index == n) {
    signal_truncated <- FALSE
    if (verbose) {
      message("  Signal length exactly matches target - no truncation needed")
    }
  }
  
  # ===========================================================================
  # Step 5: Truncate Signals
  # ===========================================================================
  
  if (verbose) message(sprintf("Truncating to %d samples...", truncation_index))
  
  trunc_x <- acc_x[1:truncation_index]
  trunc_y <- acc_y[1:truncation_index]
  trunc_z <- acc_z[1:truncation_index]
  
  # ===========================================================================
  # Step 6: Truncate Norm (pre-computed from raw signals in Step 0)
  # ===========================================================================
  
  trunc_norm <- NULL
  
  if (compute_norm) {
    if (verbose) message("Truncating pre-computed vector norm...")
    trunc_norm <- norm_full_raw[1:truncation_index]
  }
  
  # ===========================================================================
  # Step 7: Calculate Target Resampled Length
  # ===========================================================================
  
  if (is.null(target_samples)) {
    # Default: target_steps  x  samples_per_step
    target_samples <- target_steps * samples_per_step
    
    if (verbose) {
      message(sprintf("Target resampled length: %d steps  x  %d samples/step = %d samples",
                      target_steps, samples_per_step, target_samples))
    }
  } else {
    if (verbose) {
      message(sprintf("Using provided target_samples: %d", target_samples))
    }
  }
  
  # ===========================================================================
  # Step 8: Resample Signals
  # ===========================================================================
  
  if (verbose) {
    message(sprintf("Resampling: %d -> %d samples (ratio = %.4f)",
                    truncation_index, target_samples, 
                    target_samples / truncation_index))
  }
  
  # Log resampling setup to debug file (if active)
  resample_log(log_file, "Resampling: %d -> %d (ratio=%.6f, method=%s)",
               truncation_index, target_samples,
               target_samples / truncation_index, resampling_method)
  
  # Check if resample_signal function is available
  if (!exists("resample_signal", mode = "function")) {
    if (file.exists("resample_signal.R")) {
      source("resample_signal.R")
    } else {
      stop("resample_signal() function not found. Source resample_signal.R first.")
    }
  }
  
  # Resample each component
  resamp_x <- resample_signal(trunc_x, target_samples, 
                              method = resampling_method, 
                              debug = verbose,
                              log_file = log_file, log_label = "x")
  resamp_y <- resample_signal(trunc_y, target_samples, 
                              method = resampling_method,
                              debug = FALSE,
                              log_file = log_file, log_label = "y")
  resamp_z <- resample_signal(trunc_z, target_samples, 
                              method = resampling_method,
                              debug = FALSE,
                              log_file = log_file, log_label = "z")
  
  # Get resampling method that was actually used
  method_used <- attr(resamp_x, "method_used")
  was_problematic <- attr(resamp_x, "was_problematic")
  
  resamp_norm <- NULL
  if (compute_norm) {
    resamp_norm <- resample_signal(trunc_norm, target_samples,
                                   method = resampling_method,
                                   debug = FALSE,
                                   log_file = log_file, log_label = "norm")
  }
  
  if (was_problematic && verbose) {
    message("  Note: Two-stage resampling used due to ratio close to 1.0")
  }
  
  # ===========================================================================
  # Step 9: Compile Output
  # ===========================================================================
  
  # Truncated signals
  truncated <- list(
    x = trunc_x,
    y = trunc_y,
    z = trunc_z
  )
  if (compute_norm) {
    truncated$norm <- trunc_norm
  }
  
  # Resampled signals
  resampled <- list(
    x = as.numeric(resamp_x),  # Remove attributes
    y = as.numeric(resamp_y),
    z = as.numeric(resamp_z)
  )
  if (compute_norm) {
    resampled$norm <- as.numeric(resamp_norm)
  }
  
  # Metadata
  metadata <- list(
    step_freq_hz = step_freq_hz_precise,
    step_freq_bpm = step_freq_bpm,
    target_steps = target_steps,
    actual_steps = actual_steps,
    truncated_length = truncation_index,
    resampled_length = target_samples,
    samples_per_step = samples_per_step,
    sampling_freq = sampling_freq,
    resampling_ratio = target_samples / truncation_index,
    resampling_method = method_used,
    signal_truncated = signal_truncated,
    signal_too_short = signal_too_short,
    tilt_correction_applied = apply_tilt_correction
  )
  
  # Quality indicators
  quality <- list(
    
    tilt_correction_valid = tilt_correction_valid,
    warnings = warnings_list
  )
  
  if (verbose) {
    message("Processing complete.")
    message(sprintf("  Truncated: %d samples (%.1f steps)", 
                    truncation_index, actual_steps))
    message(sprintf("  Resampled: %d samples (%d samples/step)",
                    target_samples, samples_per_step))
  }
  
  # Return result
  structure(
    list(
      truncated = truncated,
      resampled = resampled,
      metadata = metadata,
      tilt_correction = tilt_diagnostics,
      quality = quality
    ),
    class = "gait_preprocessed"
  )
}


# =============================================================================
# HELPER FUNCTION: Print Method
# =============================================================================

#' Print Method for Preprocessed Gait Data
#' @export
print.gait_preprocessed <- function(x, ...) {
  cat("\n")
  cat("===============================================\n")
  cat("  Gait Signal Preprocessing Result\n")
  cat("===============================================\n\n")
  
  # Tilt correction info
  if (!is.null(x$tilt_correction)) {
    cat("Tilt Correction (Moe-Nilssen 1998):\n")
    cat(sprintf("  AP angle: %+.1f deg\n", x$tilt_correction$theta_ap_deg))
    cat(sprintf("  ML angle: %+.1f deg\n", x$tilt_correction$theta_ml_deg))
    if (x$tilt_correction$orientation_corrected) {
      cat("  [!] Orientation auto-corrected (sensor was inverted)\n")
    }
    if (!x$quality$tilt_correction_valid) {
      cat("  [!] WARNING: Quality check failed\n")
    }
    cat("\n")
  } else if (!x$metadata$tilt_correction_applied) {
    cat("Tilt Correction: Skipped (pre-corrected data)\n\n")
  }
  
  cat("Step Frequency Detection:\n")
  cat(sprintf("  %.3f Hz (%.1f steps/min)\n\n", 
              x$metadata$step_freq_hz, x$metadata$step_freq_bpm))
  
  cat("Truncated Signals (Original Rate):\n")
  cat(sprintf("  Length: %d samples (%.2f seconds at %.0f Hz)\n", 
              x$metadata$truncated_length,
              x$metadata$truncated_length / x$metadata$sampling_freq,
              x$metadata$sampling_freq))
  cat(sprintf("  Steps: %.1f ", x$metadata$actual_steps))
  if (x$metadata$signal_too_short) {
    cat(sprintf("[!] Target was %d\n", x$metadata$target_steps))
  } else {
    cat(sprintf("(target: %d) [OK]\n", x$metadata$target_steps))
  }
  cat(sprintf("  Components: %s\n\n", paste(names(x$truncated), collapse = ", ")))
  
  cat("Resampled Signals (Standardized):\n")
  cat(sprintf("  Length: %d samples\n", x$metadata$resampled_length))
  cat(sprintf("  Samples per step: %d\n", x$metadata$samples_per_step))
  cat(sprintf("  Resampling ratio: %.4f", x$metadata$resampling_ratio))
  if (x$metadata$resampling_method == "two_stage") {
    cat(" [two-stage]\n")
  } else {
    cat(" [direct]\n")
  }
  cat(sprintf("  Components: %s\n", paste(names(x$resampled), collapse = ", ")))
  
  # Warnings section
  if (length(x$quality$warnings) > 0) {
    cat("\n")
    cat("===============================================\n")
    cat("  WARNINGS\n")
    cat("===============================================\n")
    for (w in x$quality$warnings) {
      cat(sprintf("  [!] %s\n", w))
    }
  }
  
  cat("\n")
  invisible(x)
}

# =============================================================================
# HELPER FUNCTION: Summary Method
# =============================================================================

#' Summary Method for Preprocessed Gait Data
#' @export
summary.gait_preprocessed <- function(object, ...) {
  cat("Gait Preprocessing Summary\n")
  cat("--------------------------\n")
  
  # Tilt correction
  if (!is.null(object$tilt_correction)) {
    cat(sprintf("\nTilt Correction: AP=%.1f deg, ML=%.1f deg\n",
                object$tilt_correction$theta_ap_deg,
                object$tilt_correction$theta_ml_deg))
  }
  
  # Signal statistics
  cat("\nTruncated Signal Statistics:\n")
  for (comp in names(object$truncated)) {
    sig <- object$truncated[[comp]]
    cat(sprintf("  %s: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
                comp, mean(sig), sd(sig), min(sig), max(sig)))
  }
  
  cat("\nResampled Signal Statistics:\n")
  for (comp in names(object$resampled)) {
    sig <- object$resampled[[comp]]
    cat(sprintf("  %s: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
                comp, mean(sig), sd(sig), min(sig), max(sig)))
  }
  
  cat("\nKey Metadata:\n")
  cat(sprintf("  Step frequency: %.3f Hz\n", object$metadata$step_freq_hz))
  cat(sprintf("  Actual steps: %.1f\n", object$metadata$actual_steps))
  cat(sprintf("  Truncated -> Resampled: %d -> %d\n",
              object$metadata$truncated_length, object$metadata$resampled_length))
  cat(sprintf("  Tilt correction applied: %s\n", 
              ifelse(object$metadata$tilt_correction_applied, "Yes", "No")))
  
  invisible(object)
}


