#' Compute RMS and RMS Ratio for Gait Analysis
#' 
#' Computes movement intensity (RMS of vector norm) and lateral stability index
#' (RMS ratio = RMS_mediolateral / RMS_norm) from truncated lumbar accelerometer signals.
#' 
#' @description
#' Root mean square (RMS) quantifies movement intensity during walking and serves
#' as a proxy for walking speed. The RMS ratio normalizes mediolateral trunk
#' accelerations by total movement intensity, providing a speed-independent measure
#' of lateral trunk control and gait quality.
#'
#' @param az Numeric vector: mediolateral acceleration (g), TILT-CORRECTED
#' @param norm_signal Numeric vector: vector norm with gravity removed,
#'        computed from RAW signals BEFORE tilt correction:
#'        norm = sqrt(ax_raw^2 + ay_raw^2 + az_raw^2) - 1
#'
#' @return A named list with:
#'   \item{rms_norm}{RMS of vector norm (movement intensity, proxy for walking speed)}
#'   \item{rms_ratio}{RMS_ml / RMS_norm (lateral stability index)}
#'
#' @note IMPORTANT: To match the ACIER methodology, norm_signal must be computed
#'       from RAW (non-tilt-corrected) signals, while az should be tilt-corrected.
#'       This is because the norm represents total movement intensity (orientation-
#'       independent), while ML acceleration requires tilt correction for accuracy.
#'
#' @references
#' RMS methodology:
#'   Moe-Nilssen, R. (1998). A new method for evaluating motor control in gait
#'   under real-life environmental conditions. Part 2: Gait analysis.
#'   Clinical Biomechanics, 13(4-5), 328-335.
#'   https://doi.org/10.1016/S0268-0033(98)00090-4
#'
#' RMS ratio (mediolateral / vector norm):
#'   Sekine, M., Tamura, T., Yoshida, M., Suda, Y., Kimura, Y., Miyoshi, H.,
#'   Kijima, Y., Higashi, Y., & Fujimoto, T. (2013). A gait abnormality measure
#'   based on root mean square of trunk acceleration. Journal of NeuroEngineering
#'   and Rehabilitation, 10, 118. https://doi.org/10.1186/1743-0003-10-118
#'
#' ACIER study application:
#'   Piergiovanni, S., & Terrier, P. (2024). Validity of linear and nonlinear
#'   measures of gait variability to characterize aging gait with a single lower
#'   back accelerometer. Sensors, 24(23), 7427.
#'   https://doi.org/10.3390/s24237427
#'
#' @examples
#' # After tilt correction and signal truncation:
#' result <- compute_rms(az_corrected, norm_signal_raw)
#' cat("Movement intensity:", result$rms_norm, "g\n")
#' cat("Lateral stability index:", result$rms_ratio, "\n")
#'
#' @export
compute_rms <- function(az, norm_signal) {
  
  if (!is.numeric(az) || !is.numeric(norm_signal)) {
    stop("az and norm_signal must be numeric vectors")
  }
  
  if (length(az) != length(norm_signal)) {
    stop("az and norm_signal must have the same length")
  }
  
  if (length(az) == 0) {
    stop("Input signals are empty")
  }
  
  # RMS = sqrt(mean(x²))
  rms_norm <- sqrt(mean(norm_signal^2))
  rms_ml <- sqrt(mean(az^2))
  
  # RMS ratio: ML / norm
  rms_ratio <- rms_ml / rms_norm
  
  return(list(
    rms_norm = rms_norm,
    rms_ratio = rms_ratio
  ))
}
