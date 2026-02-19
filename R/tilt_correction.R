#' Moe-Nilssen (1998) Tilt Correction Algorithm
#' 
#' Transforms acceleration data from an arbitrarily tilted triaxial accelerometer
#' to a horizontal-vertical orthogonal coordinate system, removing the effects
#' of gravity and sensor misalignment.
#'
#' @description
#' This function implements the algorithm described in:
#' Moe-Nilssen, R. (1998). "A new method for evaluating motor control in gait 
#' under real-life environmental conditions. Part 1: The instrument."
#' Clinical Biomechanics, 13(4-5), 320-327. doi:10.1016/s0268-0033(98)00089-8
#'
#' The algorithm uses the accelerometer's capacity as an inclinometer: since a
#' piezoresistant accelerometer registers gravity as a static component, the mean
#' acceleration along each axis during a walking trial provides information about
#' the sensor's tilt angles relative to the horizontal plane.
#'
#' @details
#' ## Orientation Detection
#' 
#' The algorithm automatically detects sensor orientation by examining the mean
#' vertical acceleration:
#' - If mean(V) > 0: Sensor's vertical axis points UPWARD (paper's convention)
#' - If mean(V) < 0: Sensor's vertical axis points DOWNWARD (inverted mounting)
#'
#' When an inverted orientation is detected, the vertical signal is automatically
#' flipped before applying the tilt correction, ensuring consistent output
#' regardless of how the sensor was mounted.
#'
#' ## Algorithm Steps
#' 
#' 1. Estimate tilt angles from mean accelerations (gravity as inclinometer)
#' 2. First rotation: Correct AP tilt in sagittal plane (Equations A1-A2)
#' 3. Second rotation: Correct ML tilt in frontal plane (Equations A3-A4)
#' 4. Remove 1g gravity component from vertical axis
#'
#' ## Important Assumptions (from original paper)
#' 
#' - Constant tilt: Accelerometer maintains relatively constant orientation
#' - Sagittal plane alignment: AP axis lies in the sagittal plane
#' - Zero mean dynamic acceleration: Net velocity change ≈ 0 over walking trial
#' - Small to moderate tilts: Works best for typical walking postures
#' - No significant axial rotation: Rotation about vertical axis has minimal effect
#'
#' @param acc_ap Numeric vector. Anteroposterior acceleration in g units.
#' @param acc_v Numeric vector. Vertical acceleration in g units.
#' @param acc_ml Numeric vector. Mediolateral acceleration in g units.
#' @param auto_orient Logical. If TRUE (default), automatically detect and correct
#'   for inverted sensor orientation. If FALSE, assume paper's convention
#'   (vertical axis pointing upward).
#' @param return_diagnostics Logical. If TRUE, return additional diagnostic
#'   information including tilt angles and orientation detection results.
#'   Default is FALSE.
#' @param verbose Logical. If TRUE, print diagnostic messages. Default is FALSE
#'
#' @return If `return_diagnostics = FALSE` (default), returns a list with:
#' \describe{
#'   \item{ap}{Corrected AP acceleration in horizontal plane (mean ≈ 0)}
#'   \item{v}{Corrected vertical acceleration with gravity removed (mean ≈ 0)}
#'   \item{ml}{Corrected ML acceleration in horizontal plane (mean ≈ 0)}
#' }
#'
#' If `return_diagnostics = TRUE`, additionally includes:
#' \describe{
#'   \item{theta_ap}{Estimated AP tilt angle in radians}
#'   \item{theta_ml}{Estimated ML tilt angle in radians}
#'   \item{theta_ap_deg}{Estimated AP tilt angle in degrees}
#'   \item{theta_ml_deg}{Estimated ML tilt angle in degrees}
#'   \item{orientation_detected}{Character: "upward" or "downward"}
#'   \item{orientation_corrected}{Logical: was the signal flipped?}
#'   \item{mean_input}{Named vector of input mean accelerations}
#'   \item{mean_output}{Named vector of output mean accelerations}
#'   \item{quality_check}{Logical: TRUE if all output means are near zero}
#' }
#'
#' @examples
#' # Simulate tilted accelerometer data (10° forward tilt, 5° right tilt)
#' set.seed(42)
#' n <- 5000
#' theta_ap_true <- 10 * pi / 180  # 10 degrees
#' theta_ml_true <- 5 * pi / 180   # 5 degrees
#' 
#' # Pure gravity in sensor frame (tilted)
#' acc_ap <- sin(theta_ap_true) + rnorm(n, 0, 0.05)
#' acc_v <- cos(theta_ap_true) * cos(theta_ml_true) + rnorm(n, 0, 0.05)
#' acc_ml <- sin(theta_ml_true) + rnorm(n, 0, 0.05)
#' 
#' # Apply tilt correction
#' result <- tilt_correction(acc_ap, acc_v, acc_ml)
#' 
#' # Check that means are near zero
#' sapply(result, mean)
#'
#' # With diagnostics
#' result_diag <- tilt_correction(acc_ap, acc_v, acc_ml, return_diagnostics = TRUE)
#' print(result_diag$theta_ap_deg)  # Should be close to 10°
#' print(result_diag$theta_ml_deg)  # Should be close to 5°
#'
#' @references
#' Moe-Nilssen, R. (1998). A new method for evaluating motor control in gait 
#' under real-life environmental conditions. Part 1: The instrument.
#' Clinical Biomechanics, 13(4-5), 320-327. 
#' \doi{10.1016/s0268-0033(98)00089-8}
#'
#' @seealso
#' For the ACIER project context, see the Moe_Nilssen_1998_Algorithm_Analysis.md
#' documentation which details the sign conventions and implementation notes.
#'
#' @export
tilt_correction <- function(acc_ap, acc_v, acc_ml,
                            auto_orient = TRUE,
                            return_diagnostics = FALSE,
                            verbose = FALSE) {
  
  # ============================================================================
  # Input Validation
  # ============================================================================
  
  # Check that all inputs are numeric vectors of the same length
  if (!is.numeric(acc_ap) || !is.numeric(acc_v) || !is.numeric(acc_ml)) {
    stop("All acceleration inputs must be numeric vectors")
  }
  
  n <- length(acc_ap)
  if (length(acc_v) != n || length(acc_ml) != n) {
    stop("All acceleration vectors must have the same length")
  }
  
  if (n < 10) {
    warning("Very short signal (n < 10). Tilt angle estimates may be unreliable.")
  }
  
  # Check for NA values
  if (any(is.na(acc_ap)) || any(is.na(acc_v)) || any(is.na(acc_ml))) {
    warning("NA values detected in input. These will propagate to output.")
  }
  
  # ============================================================================
  # Store Original Means for Diagnostics
  # ============================================================================
  
  mean_input <- c(ap = mean(acc_ap, na.rm = TRUE),
                  v = mean(acc_v, na.rm = TRUE),
                  ml = mean(acc_ml, na.rm = TRUE))
  
  # ============================================================================
  # Orientation Detection and Correction
  # ============================================================================
  
  # The Moe-Nilssen algorithm assumes vertical axis points UPWARD
  # When standing still: mean(V) ≈ +1g
  # If sensor is mounted upside down: mean(V) ≈ -1g
  
  orientation_detected <- "upward"
  orientation_corrected <- FALSE
  
  if (auto_orient) {
    if (mean_input["v"] < -0.1) {
      # Sensor appears to be inverted (vertical pointing down)
      orientation_detected <- "downward"
      orientation_corrected <- TRUE
      
      # Flip the vertical axis to match paper's convention
      acc_v <- -acc_v
      
      if (verbose) {
        message(sprintf("Detected inverted sensor orientation (mean V = %.3f g). ",
                        mean_input["v"]),
                "Flipping vertical axis for correction.")
      }
    } else if (verbose) {
      message(sprintf("Sensor orientation appears normal (mean V = %.3f g).",
                      mean_input["v"]))
    }
  }
  
  # ============================================================================
  # Step 1: Estimate Tilt Angles from Mean Accelerations
  # ============================================================================
  # 
  # Key insight: Over a walking trial, the mean dynamic acceleration approaches
  # zero (cyclical movement). Therefore, the mean measured acceleration equals
  # the static gravity component:
  #   sin(θ_ap) ≈ mean(a_ap)
  #   sin(θ_ml) ≈ mean(a_ml)
  
  sin_theta_ap <- mean(acc_ap, na.rm = TRUE)
  sin_theta_ml <- mean(acc_ml, na.rm = TRUE)
  
  # Check that sin values are in valid range [-1, 1]
  # Clamp values slightly outside due to noise
  if (abs(sin_theta_ap) > 1) {
    if (abs(sin_theta_ap) > 1.1) {
      warning(sprintf("AP mean (%.3f g) exceeds 1g - check sensor calibration or data units.",
                      sin_theta_ap))
    }
    sin_theta_ap <- sign(sin_theta_ap) * min(abs(sin_theta_ap), 0.999)
  }
  
  if (abs(sin_theta_ml) > 1) {
    if (abs(sin_theta_ml) > 1.1) {
      warning(sprintf("ML mean (%.3f g) exceeds 1g - check sensor calibration or data units.",
                      sin_theta_ml))
    }
    sin_theta_ml <- sign(sin_theta_ml) * min(abs(sin_theta_ml), 0.999)
  }
  
  # Calculate tilt angles
  theta_ap <- asin(sin_theta_ap)
  theta_ml <- asin(sin_theta_ml)
  
  # Pre-compute trigonometric values for efficiency
  cos_theta_ap <- cos(theta_ap)
  sin_theta_ap <- sin(theta_ap)
  cos_theta_ml <- cos(theta_ml)
  sin_theta_ml <- sin(theta_ml)
  
  if (verbose) {
    message(sprintf("Estimated tilt angles: AP = %.2f°, ML = %.2f°",
                    theta_ap * 180 / pi,
                    theta_ml * 180 / pi))
  }
  
  # ============================================================================
  # Step 2: First Rotation - Correct AP Tilt (Sagittal Plane)
  # ============================================================================
  #
  # This rotation corrects for forward/backward tilt.
  # The rotation is applied to the AP and V components.
  #
  # Rotation matrix:
  # [a_A  ]   [cos(θ_ap)  -sin(θ_ap)] [a_ap]
  # [a_v' ] = [sin(θ_ap)   cos(θ_ap)] [a_v ]
  #
  # Equation A1: a_A = a_ap · cos(θ_ap) - a_v · sin(θ_ap)
  # Equation A2: a_v' = a_ap · sin(θ_ap) + a_v · cos(θ_ap)
  
  acc_A <- acc_ap * cos_theta_ap - acc_v * sin_theta_ap
  acc_v_provisional <- acc_ap * sin_theta_ap + acc_v * cos_theta_ap
  
  # ============================================================================
  # Step 3: Second Rotation - Correct ML Tilt (Frontal Plane)
  # ============================================================================
  #
  # This rotation corrects for left/right tilt.
  # The rotation is applied to the ML and provisional V components.
  #
  # Rotation matrix:
  # [a_M]   [cos(θ_ml)  -sin(θ_ml)] [a_ml ]
  # [a_V] = [sin(θ_ml)   cos(θ_ml)] [a_v' ] - [0]
  #                                           [1]
  #
  # Equation A3: a_M = a_ml · cos(θ_ml) - a_v' · sin(θ_ml)
  # Equation A4: a_V = a_ml · sin(θ_ml) + a_v' · cos(θ_ml) - 1
  #
  # Note: The -1 removes the 1g gravity component from the final vertical
  
  acc_M <- acc_ml * cos_theta_ml - acc_v_provisional * sin_theta_ml
  acc_V <- acc_ml * sin_theta_ml + acc_v_provisional * cos_theta_ml - 1
  
  # ============================================================================
  # Prepare Output
  # ============================================================================
  
  result <- list(
    ap = acc_A,
    v = acc_V,
    ml = acc_M
  )
  
  if (return_diagnostics) {
    mean_output <- c(ap = mean(acc_A, na.rm = TRUE),
                     v = mean(acc_V, na.rm = TRUE),
                     ml = mean(acc_M, na.rm = TRUE))
    
    # Quality check: output means should be near zero after correction.
    #
    # Per-axis tolerances (physically justified):
    #   AP/ML: 0.05g - the rotation equations directly cancel the gravity
    #          component on these axes; residuals > 0.05g indicate real problems
    #          (sensor malfunction, non-walking data, wrong units).
    #   V:     0.10g - depends on the zero-mean-dynamic-acceleration assumption
    #          (Moe-Nilssen 1998). Short signals, an incomplete final stride, or
    #          mildly asymmetric gait can produce residuals up to 0.10g without
    #          any actual problem. 
    quality_check <- (abs(mean_output["ap"]) < 0.05) &&
      (abs(mean_output["v"])  < 0.10) &&
      (abs(mean_output["ml"]) < 0.05)
    
    if (!quality_check && verbose) {
      message(sprintf(
        "Tilt correction quality check: AP=%.4fg  V=%.4fg  ML=%.4fg",
        mean_output["ap"], mean_output["v"], mean_output["ml"]
      ))
      if (abs(mean_output["ap"]) >= 0.05 || abs(mean_output["ml"]) >= 0.05) {
        message("  [!] AP or ML residual > 0.05g - check sensor alignment or data units")
      }
      if (abs(mean_output["v"]) >= 0.10) {
        message("  [!] V residual > 0.10g - signal may be too short or gait asymmetric")
      }
    }
    
    result$theta_ap <- theta_ap
    result$theta_ml <- theta_ml
    result$theta_ap_deg <- theta_ap * 180 / pi
    result$theta_ml_deg <- theta_ml * 180 / pi
    result$orientation_detected <- orientation_detected
    result$orientation_corrected <- orientation_corrected
    result$mean_input <- mean_input
    result$mean_output <- mean_output
    result$quality_check <- quality_check
  }
  
  return(result)
}


#' Batch Tilt Correction for Multiple Signals
#'
#' Apply tilt correction to a data frame or matrix containing triaxial
#' accelerometer data.
#'
#' @param data A data frame or matrix with acceleration data.
#' @param col_ap Column name or index for AP acceleration.
#' @param col_v Column name or index for vertical acceleration.
#' @param col_ml Column name or index for ML acceleration.
#' @param ... Additional arguments passed to `tilt_correction()`.
#'
#' @return A data frame with corrected accelerations in columns named
#'   `ap_corrected`, `v_corrected`, and `ml_corrected`.
#'
#' @examples
#' # Create sample data frame
#' df <- data.frame(
#'   time = seq(0, 10, length.out = 2560),
#'   ap = sin(0.17) + rnorm(2560, 0, 0.05),
#'   v = cos(0.17) * cos(0.09) + rnorm(2560, 0, 0.05),
#'   ml = sin(0.09) + rnorm(2560, 0, 0.05)
#' )
#'
#' # Apply correction
#' df_corrected <- tilt_correction_batch(df, "ap", "v", "ml")
#'
#' @export
tilt_correction_batch <- function(data, col_ap, col_v, col_ml, ...) {
  
  # Extract acceleration vectors
  if (is.data.frame(data)) {
    acc_ap <- data[[col_ap]]
    acc_v <- data[[col_v]]
    acc_ml <- data[[col_ml]]
  } else if (is.matrix(data)) {
    acc_ap <- data[, col_ap]
    acc_v <- data[, col_v]
    acc_ml <- data[, col_ml]
  } else {
    stop("'data' must be a data.frame or matrix")
  }
  
  # Apply tilt correction
  corrected <- tilt_correction(acc_ap, acc_v, acc_ml, ...)
  
  # Create output data frame
  result <- data.frame(
    ap_corrected = corrected$ap,
    v_corrected = corrected$v,
    ml_corrected = corrected$ml
  )
  
  return(result)
}


#' Validate Tilt Correction Results
#'
#' Check that tilt correction was successful by verifying that the output
#' signals have means near zero and that tilt angles are within reasonable
#' ranges.
#'
#' @param corrected_data Output from `tilt_correction()` with
#'   `return_diagnostics = TRUE`.
#' @param mean_tolerance Maximum acceptable mean acceleration (in g) for
#'   quality check. Default is 0.01g.
#' @param max_tilt_deg Maximum acceptable tilt angle (in degrees) for
#'   warning. Default is 30°.
#'
#' @return A list with validation results:
#' \describe{
#'   \item{valid}{Logical: overall validation passed?}
#'   \item{mean_check}{Logical: all means within tolerance?}
#'   \item{tilt_check}{Logical: tilt angles within reasonable range?}
#'   \item{messages}{Character vector of warnings/messages}
#' }
#'
#' @export
validate_tilt_correction <- function(corrected_data,
                                     mean_tolerance_ap = 0.05,
                                     mean_tolerance_v  = 0.10,
                                     mean_tolerance_ml = 0.05,
                                     max_tilt_deg = 30) {
  
  if (is.null(corrected_data$mean_output)) {
    stop("Input must be from tilt_correction() with return_diagnostics = TRUE")
  }
  
  messages <- character(0)
  
  # Per-axis mean checks (tolerances match tilt_correction() quality_check)
  ap_ok <- abs(corrected_data$mean_output["ap"]) < mean_tolerance_ap
  v_ok  <- abs(corrected_data$mean_output["v"])  < mean_tolerance_v
  ml_ok <- abs(corrected_data$mean_output["ml"]) < mean_tolerance_ml
  mean_check <- ap_ok && v_ok && ml_ok
  
  if (!ap_ok) {
    messages <- c(messages, sprintf(
      "AP mean too large: %.4fg (tolerance: %.4fg) - check sensor alignment",
      corrected_data$mean_output["ap"], mean_tolerance_ap))
  }
  if (!v_ok) {
    messages <- c(messages, sprintf(
      "V mean too large: %.4fg (tolerance: %.4fg) - signal too short or asymmetric gait",
      corrected_data$mean_output["v"], mean_tolerance_v))
  }
  if (!ml_ok) {
    messages <- c(messages, sprintf(
      "ML mean too large: %.4fg (tolerance: %.4fg) - check sensor alignment",
      corrected_data$mean_output["ml"], mean_tolerance_ml))
  }
  
  # Check tilt angles
  tilt_check <- abs(corrected_data$theta_ap_deg) <= max_tilt_deg &&
    abs(corrected_data$theta_ml_deg) <= max_tilt_deg
  
  if (!tilt_check) {
    messages <- c(messages,
                  sprintf("Extreme tilt detected: AP=%.1f°, ML=%.1f° (max=%.1f°)",
                          corrected_data$theta_ap_deg,
                          corrected_data$theta_ml_deg,
                          max_tilt_deg))
  }
  
  # Check orientation correction
  if (corrected_data$orientation_corrected) {
    messages <- c(messages, "Sensor orientation was automatically corrected (inverted mounting detected)")
  }
  
  list(
    valid = mean_check && tilt_check,
    mean_check = mean_check,
    tilt_check = tilt_check,
    messages = messages
  )
}


# ==============================================================================
# UNIT TESTS / EXAMPLES
# ==============================================================================

if (FALSE) {  # Don't run automatically
  
  # ---------------------------------------------------------------------------
  # Test 1: Known tilt angles
  # ---------------------------------------------------------------------------
  
  cat("\n=== Test 1: Known tilt angles ===\n")
  
  set.seed(42)
  n <- 10000
  
  # True tilt: 10° forward, 5° right
  theta_ap_true <- 10 * pi / 180
  theta_ml_true <- 5 * pi / 180
  
  # Generate "sensor" accelerations (gravity decomposed)
  # When tilted, gravity appears in AP and ML channels
  acc_ap <- sin(theta_ap_true) + rnorm(n, 0, 0.03)  # ~0.17g mean
  acc_v <- cos(theta_ap_true) * cos(theta_ml_true) + rnorm(n, 0, 0.03)  # ~0.98g mean
  acc_ml <- sin(theta_ml_true) + rnorm(n, 0, 0.03)  # ~0.09g mean
  
  result <- tilt_correction(acc_ap, acc_v, acc_ml, 
                            return_diagnostics = TRUE, 
                            verbose = TRUE)
  
  cat("\nEstimated tilt angles:\n")
  cat(sprintf("  AP: %.2f° (true: %.2f°)\n", result$theta_ap_deg, 10))
  cat(sprintf("  ML: %.2f° (true: %.2f°)\n", result$theta_ml_deg, 5))
  
  cat("\nOutput means (should be ~0):\n")
  print(result$mean_output)
  
  cat("\nQuality check:", result$quality_check, "\n")
  
  # ---------------------------------------------------------------------------
  # Test 2: Inverted sensor (vertical pointing down)
  # ---------------------------------------------------------------------------
  
  cat("\n=== Test 2: Inverted sensor orientation ===\n")
  
  # Same tilt, but vertical axis inverted
  acc_v_inverted <- -acc_v
  
  result_inv <- tilt_correction(acc_ap, acc_v_inverted, acc_ml,
                                return_diagnostics = TRUE,
                                verbose = TRUE)
  
  cat("\nOrientation detected:", result_inv$orientation_detected, "\n")
  cat("Orientation corrected:", result_inv$orientation_corrected, "\n")
  cat("\nOutput means (should still be ~0):\n")
  print(result_inv$mean_output)
  
  # ---------------------------------------------------------------------------
  # Test 3: Zero tilt (ideal mounting)
  # ---------------------------------------------------------------------------
  
  cat("\n=== Test 3: Zero tilt (ideal mounting) ===\n")
  
  acc_ap_ideal <- rnorm(n, 0, 0.05)
  acc_v_ideal <- 1 + rnorm(n, 0, 0.05)
  acc_ml_ideal <- rnorm(n, 0, 0.05)
  
  result_ideal <- tilt_correction(acc_ap_ideal, acc_v_ideal, acc_ml_ideal,
                                  return_diagnostics = TRUE,
                                  verbose = TRUE)
  
  cat("\nEstimated tilt angles (should be ~0):\n")
  cat(sprintf("  AP: %.2f°\n", result_ideal$theta_ap_deg))
  cat(sprintf("  ML: %.2f°\n", result_ideal$theta_ml_deg))
  
  # ---------------------------------------------------------------------------
  # Test 4: Validation function
  # ---------------------------------------------------------------------------
  
  cat("\n=== Test 4: Validation ===\n")
  
  validation <- validate_tilt_correction(result)
  print(validation)
}
