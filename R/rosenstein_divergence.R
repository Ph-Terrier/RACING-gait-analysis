#' =============================================================================
#' Rosenstein Divergence Curve Computation
#' =============================================================================
#' 
#' Computes divergence curves from time series using Rosenstein's algorithm
#' for estimating the largest Lyapunov exponent. This implementation focuses
#' on computing the raw divergence curve; slope fitting for LDS/ACI is performed
#' separately using domain-specific stride ranges (e.g., 0-0.5, 5-12 strides).
#'
#' Scientific Reference:
#'   Rosenstein, M. T., Collins, J. J., & De Luca, C. J. (1993).
#'   A practical method for calculating largest Lyapunov exponents from small data sets.
#'   Physica D: Nonlinear Phenomena, 65(1-2), 117-134.
#'   https://doi.org/10.1016/0167-2789(93)90009-P
#'
#' MATLAB Implementation Source:
#'   This R implementation is adapted from MATLAB code originally developed by:
#'   mirwais (2013). "Largest Lyapunov Exponent with Rosenstein's Algorithm"
#'   MATLAB Central File Exchange.
#'   https://www.mathworks.com/matlabcentral/fileexchange/38424
#'   (Retrieved February 2025)
#'
#' Alternative Reference Implementation:
#'   TISEAN 3.0 - Nonlinear Time Series Analysis
#'   Hegger, R., Kantz, H., & Schreiber, T. (1999).
#'   https://www.pks.mpg.de/tisean/
#'
#' Key Algorithm Steps (Rosenstein et al. 1993):
#'   1. Phase space reconstruction via Method of Delays (Takens' theorem)
#'   2. Nearest neighbor search with temporal exclusion (Theiler window)
#'   3. Track divergence of nearby trajectories over time
#'   4. Average logarithmic divergence at each time step (Eq. 13 in paper)
#'
#' R Implementation Enhancements:
#'   - Fast KD-tree nearest neighbor search via RANN package (optional)
#'   - Vectorized divergence computation for performance
#'   - Graceful fallback to pure R if RANN unavailable
#'   - Robust mean period estimation with multiple fallback strategies
#'
#' ACIER Study Application:
#'   In the RACING/ACIER gait analysis pipeline, this function computes
#'   divergence curves up to ~1800 samples (12 strides). Separate fitting
#'   functions then extract:
#'   - LDS (Local Dynamic Stability): slope over 0-75 samples (0-0.5 strides)
#'   - ACI (Attractor Complexity Index): slope over 750-1800 samples (5-12 strides)
#'
#' Author: Philippe Terrier, RACING Project (ORD 2025), HE-Arc 
#' =============================================================================

#' Phase Space Reconstruction using Method of Delays (Takens' Embedding)
#' 
#' Implements Takens' embedding theorem to reconstruct the phase space
#' from a univariate time series. Each row of the output matrix represents
#' a delay vector: [x(t), x(t+tau), x(t+2*tau), ..., x(t+(m-1)*tau)].
#'
#' @references
#' Takens, F. (1981). Detecting strange attractors in turbulence.
#' In Dynamical Systems and Turbulence, Lecture Notes in Mathematics, 898, 366-381.
#' Springer, Berlin. 
#' @param x Numeric vector - the time series
#' @param m Integer - embedding dimension
#' @param tau Integer - time delay (lag)
#' @return Matrix of dimension M x m where M = length(x) - (m-1)*tau
embed_time_series <- function(x, m, tau) {
  N <- length(x)
  M <- N - (m - 1) * tau
  
  if (M <= 0) {
    stop("Time series too short for given embedding parameters")
  }
  
  indices <- outer(1:M, seq(0, (m - 1) * tau, by = tau), "+")
  matrix(x[indices], nrow = M, ncol = m)
}


#' Find Nearest Neighbors with Temporal Exclusion (Theiler Window)
#' 
#' For each point in the reconstructed phase space, finds the nearest neighbor
#' excluding points within a temporal window of +/- mean_period. This "Theiler window"
#' prevents autocorrelation effects from dominating the nearest neighbor selection.
#' 
#' Performance: Uses RANN package KD-tree (O(M log M)) when available,
#' otherwise falls back to vectorized R implementation (O(M^2)).
#'
#' @references
#' Theiler, J. (1986). Spurious dimension from correlation algorithms applied
#' to limited time-series data. Physical Review A, 34(3), 2427-2432.
#'
#' RANN: Arya, S., Mount, D. M., et al. (2019). RANN: Fast Nearest Neighbour Search.
#' R package. https://CRAN.R-project.org/package=RANN
#' 
#' @param Y Matrix - embedded phase space (M x m)
#' @param mean_period Integer - minimum temporal separation required
#' @return Integer vector of length M with indices of nearest neighbors
find_nearest_neighbors <- function(Y, mean_period, verbose = FALSE) {
  M <- nrow(Y)
  
  # Try fast KD-tree search if RANN available (significantly faster for large M)
  use_rann <- requireNamespace("RANN", quietly = TRUE)
  
  if (use_rann && verbose) {
    message("Using RANN KD-tree for fast neighbor search")
  }
  
  if (use_rann)  {
    k <- min(2 * mean_period + 10, M - 1)
    nn_result <- RANN::nn2(Y, Y, k = k + 1)
    
    nearest_pos <- integer(M)
    for (i in seq_len(M)) {
      candidates <- nn_result$nn.idx[i, -1]
      valid <- which(abs(candidates - i) > mean_period)
      
      if (length(valid) > 0) {
        nearest_pos[i] <- candidates[valid[1]]
      } else {
        # Fallback for this point
        distances <- sqrt(rowSums((Y - matrix(Y[i, ], M, ncol(Y), byrow = TRUE))^2))
        distances[max(1, i - mean_period):min(M, i + mean_period)] <- Inf
        nearest_pos[i] <- which.min(distances)
      }
    }
    return(nearest_pos)
  }
  
  # Vectorized fallback
  nearest_pos <- integer(M)
  for (i in seq_len(M)) {
    diff_mat <- Y - matrix(Y[i, ], M, ncol(Y), byrow = TRUE)
    distances <- sqrt(rowSums(diff_mat^2))
    distances[i] <- Inf
    distances[max(1, i - mean_period):min(M, i + mean_period)] <- Inf
    nearest_pos[i] <- which.min(distances)
  }
  nearest_pos
}


#' Compute Divergence Curve via Logarithmic Distance Tracking
#' 
#' Implements Equation 13 from Rosenstein et al. (1993):
#' d(k) = <ln(d_j(k))> where d_j(k) = ||Y(j+k) - Y(j*+k)||
#' 
#' For each pair of initial nearest neighbors (j, j*), tracks how their
#' Euclidean distance evolves over k time steps. The average logarithmic
#' divergence at each step k provides the raw divergence curve.
#'
#' @details
#' Zero distances are excluded to avoid -Inf in logarithm. Pairs where
#' either trajectory would exceed signal length are excluded from averaging.
#' The slope of d(k) vs k (in the appropriate range) gives the largest
#' Lyapunov exponent.
#' @param Y Matrix - embedded phase space (M x m)
#' @param nearest_pos Integer vector - indices of nearest neighbors
#' @param max_iter Integer - maximum number of time steps to track divergence
#' @return Numeric vector of average log-divergence at each time step
compute_divergence <- function(Y, nearest_pos, max_iter) {
  M <- nrow(Y)
  divergence <- numeric(max_iter)
  
  for (k in seq_len(max_iter)) {
    max_ind <- M - k
    valid <- which(seq_len(M) <= max_ind & nearest_pos <= max_ind)
    
    if (length(valid) == 0) {
      divergence[k] <- NA
      next
    }
    
    j_evolved <- valid + k
    neighbor_evolved <- nearest_pos[valid] + k
    
    diff_mat <- Y[j_evolved, , drop = FALSE] - Y[neighbor_evolved, , drop = FALSE]
    dist_k <- sqrt(rowSums(diff_mat^2))
    
    nonzero <- dist_k > 0
    divergence[k] <- if (sum(nonzero) > 0) mean(log(dist_k[nonzero])) else NA
  }
  
  divergence
}


#' Estimate Mean Period from Time Series
#' 
#' Estimates the dominant periodicity for use as the Theiler window.
#' Strategy: (1) first zero crossing of ACF, (2) first local minimum,
#' (3) fallback to N/20. This provides a conservative estimate suitable
#' for most gait signals where stride period is the dominant frequency.
#' 
#' @param x Numeric vector - the time series
#' @return Integer - estimated mean period in samples
#' @note For gait analysis, this typically returns ~150 samples (one stride)
#'       at 256 Hz sampling frequency.
estimate_mean_period <- function(x) {
  N <- length(x)
  max_lag <- min(N - 1, 500)
  
  acf_vals <- acf(x, lag.max = max_lag, plot = FALSE)$acf[, 1, 1]
  
  # Find first zero crossing
  zero_cross <- which(diff(sign(acf_vals)) != 0)[1]
  if (!is.na(zero_cross)) return(zero_cross)
  
  # Fallback: first minimum
  local_min <- which(diff(sign(diff(acf_vals))) > 0) + 1
  if (length(local_min) > 0) return(local_min[1])
  
  # Default fallback
  round(N / 20)
}


#' Compute Rosenstein Divergence Curve
#' 
#' Main function for computing the divergence curve from a time series.
#' The output can be used for fitting log-divergence slopes at various
#' stride ranges typical in gait analysis.
#' 
#' @param x Numeric vector - input time series
#' @param m Integer - embedding dimension (default: 5)
#' @param tau Integer - time delay in samples (default: 10)
#' @param mean_period Integer - temporal exclusion window in samples.
#'                    If NULL, estimated from autocorrelation.
#' @param max_iter Integer - maximum iterations for divergence tracking.
#'                 For gait analysis, set based on stride duration and
#'                 desired analysis range (e.g., 10 strides).
#' @param verbose Logical - print progress messages (default: FALSE)
#' @return List containing:
#'   - divergence: Numeric vector of average log-divergence at each time step
#'   - time_steps: Integer vector from 1 to length(divergence)
#'   - params: List of parameters used (m, tau, mean_period, N, M)
#' 
#' @examples
#' # Basic usage
#' result <- rosenstein_divergence(x, m = 5, tau = 10, max_iter = 1800)
#' 
#' # For gait analysis with known stride duration (e.g., 100 samples/stride)
#' # To analyze up to 10 strides:
#' result <- rosenstein_divergence(x, m = 5, tau = 10, max_iter = 1000)
#' 
#' # Divergence curve for external fitting
#' plot(result$time_steps, result$divergence, type = "l",
#'      xlab = "Time steps (samples)", ylab = "Mean log divergence")
#'      
rosenstein_divergence <- function(x, m = 5, tau = 10, mean_period = NULL,
                                   max_iter = 1800, verbose = FALSE) {
  x <- as.numeric(x)
  N <- length(x)
  
  # Estimate mean period if not provided
  if (is.null(mean_period)) {
    mean_period <- estimate_mean_period(x)
    if (verbose) message(sprintf("Estimated mean period: %d samples", mean_period))
  }
  
  # Phase space reconstruction
  Y <- embed_time_series(x, m, tau)
  M <- nrow(Y)
  max_iter <- min(max_iter, M - 1)
  
  if (verbose) message(sprintf("Phase space: %d vectors of dimension %d", M, m))
  
  # Find nearest neighbors
  if (verbose) message("Finding nearest neighbors...")
  nearest_pos <- find_nearest_neighbors(Y, mean_period, verbose) 
  
  # Compute divergence curve
  if (verbose) message("Computing divergence curve...")
  divergence <- compute_divergence(Y, nearest_pos, max_iter)
  
  list(
    divergence = divergence,
    time_steps = seq_len(length(divergence)),
    params = list(
      m = m,
      tau = tau,
      mean_period = mean_period,
      max_iter = max_iter,
      N = N,
      M = M
    )
  )
}



