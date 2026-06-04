# =============================================================================
# AMI - Average Mutual Information (MATLAB-Faithful Implementation)
# =============================================================================
# 
# R port of MATLAB AMI functions for the RACING/ORD2025 project.
# Original MATLAB code by Durga Lal Shrestha (2005).
#
# Purpose: Determine optimal time delay (tau) for phase space reconstruction
# before computing Lyapunov exponents via Rosenstein's algorithm.
#
# The first minimum of the AMI curve indicates the optimal embedding delay,
# which is more robust than autocorrelation zero-crossing methods.
#
# Original MATLAB code adapted from:
#   Shrestha, D. L. (2005). Average Mutual Information [MATLAB code].
#   Files: ami.m, prob.m, probxy.m, rhist.m
#   Email: durgals@hotmail.com
#   Version: 1.1.0 (June 27 - July 1, 2005)
#   Copyright 2004-2005 by Durga Lal Shrestha
#
# Reference: Fraser, A.M. & Swinney, H.L. (1986). Independent coordinates
#            for strange attractors from mutual information.
#            Physical Review A, 33(2), 1134-1140.
#            https://doi.org/10.1103/PhysRevA.33.1134
#
# Typical ACIER pipeline usage: ami(signal, n_bins = 64, n_lags = 100)
#
# IMPORTANT - Histogram binning note (2026-03):
# The original MATLAB code has an internal inconsistency:
#   - prob.m  -> rhist.m -> hist(y,M): CENTER-based, bin_width = range/(M-1)
#   - probxy.m -> computeEdge()+histc(): EDGE-based, bin_width = range/M
# This R port deliberately replicates this inconsistency to ensure AMI values
# match MATLAB output exactly. rhist() uses center-based (range/(M-1)) for
# 1D marginals, while prob_2d() uses edge-based (range/M) for joint probs.
# =============================================================================


# -----------------------------------------------------------------------------
# Helper: Relative histogram (MATLAB hist(y,M) compatible)
# -----------------------------------------------------------------------------
#' Compute relative frequency histogram replicating MATLAB hist(y, M)
#'
#' MATLAB hist(y, M) creates M bin CENTERS via linspace(min, max, M),
#' giving bin_width = range/(M-1). This differs from the edge-based
#' approach (range/M) used in probxy.m for joint probabilities.
#' We deliberately replicate this behavior for MATLAB compatibility
#' in the 1D marginal probability calculation (prob_1d).
#'
#' @param y Numeric vector
#' @param n_bins Number of bins (default 10)
#' @return List with components:
#'   - freq: Relative frequencies for each bin (length = n_bins)
#'   - centers: Bin centers (matching MATLAB linspace output)
#'   - n_bins: Number of bins used
rhist <- function(y, n_bins = 10) {
  y <- as.numeric(y)
  n <- length(y)
  n_bins <- as.integer(n_bins)
  
  min_y <- min(y)
  max_y <- max(y)
  
  # Handle edge case where all values are the same
  if (min_y == max_y) {
    freq <- rep(0, n_bins)
    freq[ceiling(n_bins / 2)] <- 1
    centers <- rep(min_y, n_bins)
    return(list(freq = freq, centers = centers, n_bins = n_bins))
  }
  
  # --- MATLAB hist(y, M) CENTER-BASED binning ---
  # Centers: linspace(min, max, M) -> bin_width = range / (M-1)
  # NOT range/M as in edge-based approaches
  centers <- seq(min_y, max_y, length.out = n_bins)
  
  # Bin edges at midpoints between consecutive centers
  # First bin extends from -Inf, last bin extends to +Inf
  midpoints <- (centers[-n_bins] + centers[-1]) / 2
  edges <- c(-Inf, midpoints, Inf)
  
  # Assign values to bins using findInterval
  # With -Inf/Inf boundaries and rightmost.closed=TRUE,
  # this replicates MATLAB hist nearest-center assignment
  bin_idx <- findInterval(y, edges, rightmost.closed = TRUE)
  
  # Safety clamp (should not be needed with -Inf/Inf boundaries)
  bin_idx[bin_idx < 1] <- 1
  bin_idx[bin_idx > n_bins] <- n_bins
  
  # Count occurrences in each bin
  counts <- tabulate(bin_idx, nbins = n_bins)
  
  list(
    freq = counts / n,
    centers = centers,
    n_bins = n_bins
  )
}


# -----------------------------------------------------------------------------
# Helper: Check for zero-probability bins
# -----------------------------------------------------------------------------
#' Check if any bin has zero probability
#' 
#' @param y Numeric vector
#' @param n_bins Number of bins
#' @return Logical: TRUE if any bin is empty
has_zero_bins <- function(y, n_bins) {
  result <- rhist(y, n_bins)
  any(result$freq == 0)
}


# -----------------------------------------------------------------------------
# 1D Probability Distribution (MATLAB prob equivalent)
# -----------------------------------------------------------------------------
#' Compute 1D probability distribution with adaptive binning
#' 
#' Computes probability distribution, automatically reducing bin count
#' if any bins would be empty. Uses bisection algorithm matching
#' MATLAB prob() behavior exactly.
#'
#' @param y Numeric vector
#' @param max_bins Maximum number of bins (default 10)
#' @return List with components:
#'   - prob: Probability distribution (vector of length n_bins)
#'   - n_bins: Actual number of bins used
prob_1d <- function(y, max_bins = 10) {
  y <- as.numeric(y)
  max_bins <- floor(max_bins)
  
  # Check first iteration - if no empty bins, accept immediately
  # This matches MATLAB prob.m line 52-54
  rh <- rhist(y, max_bins)
  if (!any(rh$freq == 0)) {
    return(list(prob = rh$freq, n_bins = rh$n_bins))
  }
  
  # Adaptive bin search using bisection (matches MATLAB exactly)
  pre_bin <- 0
  is_nonzero_bin <- FALSE
  current_bins <- max_bins
  zero_bin <- max_bins
  nonzero_bin <- max_bins
  
  while (pre_bin != current_bins) {
    has_zero <- has_zero_bins(y, current_bins)
    
    if (!has_zero) {
      # No zero distribution: record and try more bins
      tmp_bin <- current_bins
      nonzero_bin <- current_bins
      current_bins <- floor((zero_bin + nonzero_bin) / 2)
      pre_bin <- tmp_bin
      is_nonzero_bin <- TRUE
    } else {
      # Zero distribution found: reduce bins
      if (!is_nonzero_bin) {
        # Previous also had zeros: halve
        pre_bin <- current_bins
        zero_bin <- current_bins
        current_bins <- floor(current_bins / 2)
      } else {
        # Previous was okay: bisect
        tmp_bin <- current_bins
        zero_bin <- current_bins
        current_bins <- floor((zero_bin + nonzero_bin) / 2)
        pre_bin <- tmp_bin
      }
    }
    
    # Safety: prevent infinite loop
    if (current_bins < 2) {
      current_bins <- 2
      break
    }
  }
  
  # Compute final histogram
  result <- rhist(y, current_bins)
  list(prob = result$freq, n_bins = result$n_bins)
}


# -----------------------------------------------------------------------------
# 2D Joint Probability Distribution (MATLAB probxy equivalent)
# -----------------------------------------------------------------------------
#' Compute 2D joint probability distribution
#' 
#' Uses MATLAB-style edge computation with -Inf/Inf boundaries and
#' left-closed interval binning.
#' 
#' @param x Numeric vector (first variable)
#' @param y Numeric vector (second variable, same length as x)
#' @param n_bins_x Number of bins for x
#' @param n_bins_y Number of bins for y
#' @return Matrix of joint probabilities (n_bins_x x n_bins_y)
prob_2d <- function(x, y, n_bins_x = 10, n_bins_y = 10) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)
  n_bins_x <- as.integer(n_bins_x)
  n_bins_y <- as.integer(n_bins_y)
  
  if (length(y) != n) {
    stop("x and y must have the same length")
  }
  
  # Compute data ranges
  min_x <- min(x); max_x <- max(x)
  min_y <- min(y); max_y <- max(y)
  
  # Handle edge cases
  if (min_x == max_x) max_x <- min_x + 1
  if (min_y == max_y) max_y <- min_y + 1
  
  bw_x <- (max_x - min_x) / n_bins_x
  bw_y <- (max_y - min_y) / n_bins_y
  
  # MATLAB-style edges: [-Inf, edge(2:n-1), Inf]
  # This matches computeEdge() in probxy.m
  all_edges_x <- min_x + bw_x * (0:n_bins_x)
  edges_x <- c(-Inf, all_edges_x[2:n_bins_x], Inf)
  
  all_edges_y <- min_y + bw_y * (0:n_bins_y)
  edges_y <- c(-Inf, all_edges_y[2:n_bins_y], Inf)
  
  # Initialize count matrix
  nn <- matrix(0, nrow = n_bins_x, ncol = n_bins_y)
  
  # MATLAB probxy loop: X >= edgeX(i) & X < edgeX(i+1)
  for (i in 1:n_bins_x) {
    # Find X values in this bin (left-closed, right-open)
    ind_x <- which(x >= edges_x[i] & x < edges_x[i + 1])
    
    if (length(ind_x) > 0) {
      y_found <- y[ind_x]
      
      # MATLAB histc behavior for Y
      # histc: counts in [edge(k), edge(k+1)) except last bin includes edge(end)
      for (j in 1:n_bins_y) {
        if (j < n_bins_y) {
          count_j <- sum(y_found >= edges_y[j] & y_found < edges_y[j + 1])
        } else {
          # Last bin includes the right boundary (matches histc behavior)
          count_j <- sum(y_found >= edges_y[j])
        }
        nn[i, j] <- count_j
      }
    }
  }
  
  # Return as probability matrix
  nn / n
}


# -----------------------------------------------------------------------------
# Mutual Information at Single Lag
# -----------------------------------------------------------------------------
#' Compute mutual information between signal and lagged version
#' 
#' @param x Numeric vector (original signal portion)
#' @param y Numeric vector (lagged signal portion, same length)
#' @param n_bins_x Number of bins for x
#' @param n_bins_y Number of bins for y
#' @param px Probability distribution of x
#' @param py Probability distribution of y
#' @return Scalar mutual information value (bits)
mutual_info <- function(x, y, n_bins_x, n_bins_y, px, py) {
  # Joint probability
  pxy <- prob_2d(x, y, n_bins_x, n_bins_y)
  
  # Compute MI: sum p(x,y) * log2(p(x,y) / (p(x) * p(y)))
  mi <- 0
  for (i in seq_len(n_bins_x)) {
    for (j in seq_len(n_bins_y)) {
      if (pxy[i, j] > 0) {
        mi <- mi + pxy[i, j] * log2(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }
  
  mi
}


# -----------------------------------------------------------------------------
# Main AMI Function
# -----------------------------------------------------------------------------
#' Compute Average Mutual Information and Correlation for Lagged Time Series
#' 
#' Computes AMI and correlation for lags 0 to n_lags. The first minimum
#' of the AMI curve indicates the optimal embedding delay (tau) for
#' phase space reconstruction.
#'
#' @param signal Numeric vector (univariate time series) or matrix 
#'   (bivariate: columns are x and y)
#' @param n_bins Number of bins for probability estimation (default 64)
#' @param n_lags Number of lags to compute (default 100)
#' @param plot Logical: create diagnostic plot? (default FALSE)
#' @return List with components:
#'   - ami: Vector of AMI values for lags 0:n_lags (length n_lags + 1)
#'   - correlation: Vector of correlation values for lags 0:n_lags
#'   - lags: Vector of lag values (0:n_lags)
#'   - optimal_lag: First minimum of AMI curve (optimal tau)
#' @examples
#' # Typical usage for gait analysis:
#' signal <- rnorm(18750)  # Standardized 250 steps
#' result <- ami(signal, n_bins = 64, n_lags = 100)
#' tau <- result$optimal_lag
#' 
#' @export
ami <- function(signal, n_bins = 64, n_lags = 100, plot = FALSE) {
  
  # --- Input validation ---
  if (is.matrix(signal) || is.data.frame(signal)) {
    signal <- as.matrix(signal)
    if (ncol(signal) > nrow(signal)) {
      signal <- t(signal)
    }
    if (ncol(signal) == 2) {
      # Bivariate case
      x <- signal[, 1]
      y <- signal[, 2]
    } else if (ncol(signal) == 1) {
      # Univariate as matrix
      x <- signal[, 1]
      y <- x
    } else {
      stop("Signal must be univariate or bivariate (1 or 2 columns)")
    }
  } else {
    # Univariate vector
    x <- as.numeric(signal)
    y <- x
  }
  
  n <- length(x)
  n_bins <- floor(n_bins)
  n_lags <- floor(n_lags)
  
  if (n_lags < 0) {
    stop("n_lags must be non-negative")
  }
  if (n_lags >= n) {
    stop("n_lags must be less than signal length")
  }
  if (n_bins < 2) {
    stop("n_bins must be at least 2")
  }
  
  # --- Compute AMI and correlation for each lag ---
  ami_values <- numeric(n_lags + 1)
  corr_values <- numeric(n_lags + 1)
  
  for (lag in 0:n_lags) {
    # Create lagged versions (matches MATLAB indexing)
    x_lag <- x[1:(n - lag)]
    y_lag <- y[(lag + 1):n]
    
    # Compute 1D probabilities with adaptive binning
    px_result <- prob_1d(x_lag, n_bins)
    py_result <- prob_1d(y_lag, n_bins)
    
    # Compute mutual information
    ami_values[lag + 1] <- mutual_info(
      x_lag, y_lag,
      px_result$n_bins, py_result$n_bins,
      px_result$prob, py_result$prob
    )
    
    # Compute correlation
    corr_values[lag + 1] <- cor(x_lag, y_lag)
  }
  
  # --- Find first minimum of AMI curve ---
  optimal_lag <- find_first_minimum(ami_values)
  
  # --- Optional plot ---
  if (plot && n_lags > 1) {
    ami_plot(0:n_lags, ami_values, corr_values)
  }
  
  # --- Return results ---
  list(
    ami = ami_values,
    correlation = corr_values,
    lags = 0:n_lags,
    optimal_lag = optimal_lag
  )
}


# -----------------------------------------------------------------------------
# Find First Minimum (Optimal Delay Detection)
# -----------------------------------------------------------------------------
#' Find first local minimum in AMI curve
#' 
#' The first minimum of the AMI curve indicates the optimal embedding
#' delay for phase space reconstruction.
#'
#' @param ami_values Numeric vector of AMI values
#' @return Integer: MATLAB-compatible tau (1-based position in AMI curve).
#'   This matches MATLAB's find(diff(ami)>0)(1) convention used in ACIER.
#'   Note: the true 0-based lag of the minimum = return_value - 1.
#' @export
find_first_minimum <- function(ami_values) {
  n <- length(ami_values)
  
  if (n < 3) {
    return(NA_integer_)
  }
  
  # Handle NA values
  if (any(is.na(ami_values))) {
    warning("NA values in AMI curve, results may be unreliable")
    ami_values[is.na(ami_values)] <- Inf
  }
  
  # Find first local minimum
  for (i in 2:(n - 1)) {
    if (!is.na(ami_values[i]) && !is.na(ami_values[i - 1]) && !is.na(ami_values[i + 1])) {
      if (ami_values[i] < ami_values[i - 1] && 
          ami_values[i] < ami_values[i + 1]) {
        return(i)  
      }
    }
  }
  
  # No minimum found: return 1-based position where AMI drops to 1/e of initial
  # (consistent with main return convention = MATLAB find(diff>0)(1))
  threshold <- ami_values[1] / exp(1)
  below_threshold <- which(ami_values < threshold & !is.na(ami_values))
  
  if (length(below_threshold) > 0) {
    return(below_threshold[1])  # 1-based position (matches main return convention)
  }
  
  # Fallback: return NA
  NA_integer_
}


# -----------------------------------------------------------------------------
# Diagnostic Plot
# -----------------------------------------------------------------------------
#' Plot AMI and correlation curves
#' 
#' @param lags Vector of lag values
#' @param ami_values Vector of AMI values
#' @param corr_values Vector of correlation values
ami_plot <- function(lags, ami_values, corr_values) {
  # Save and restore graphics parameters (RStudio-safe margins)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mar = c(4, 4, 3, 4))
  
  # Plot correlation on left axis
  plot(lags, corr_values, type = "l", lwd = 2, col = "blue",
       xlab = "Time Lag", ylab = "Correlation",
       main = "Average Mutual Information and Correlation",
       ylim = c(min(corr_values), 1))
  
  # Add AMI on right axis
  par(new = TRUE)
  plot(lags, ami_values, type = "l", lwd = 2, col = "red",
       axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext("AMI (bits)", side = 4, line = 2.5)
  
  # Mark first minimum
  opt_lag <- find_first_minimum(ami_values)
  if (!is.na(opt_lag)) {
    # opt_lag is 1-based index in ami_values (= MATLAB convention for tau)
    # Actual lag on x-axis = opt_lag - 1 (since ami_values[1] = lag 0)
    lag_value <- opt_lag - 1
    points(lag_value, ami_values[opt_lag], pch = 19, col = "red", cex = 1.5)
    abline(v = lag_value, lty = 2, col = "gray")
    text(lag_value, ami_values[opt_lag], 
         labels = sprintf("tau = %d (lag %d)", opt_lag, lag_value),
         pos = 4, col = "red")
  }
  
  legend("topright", legend = c("Correlation", "AMI"),
         col = c("blue", "red"), lwd = 2, bty = "n")
}


# -----------------------------------------------------------------------------
# Convenience Function: Estimate Optimal Time Delay
# -----------------------------------------------------------------------------
#' Estimate optimal time delay for phase space reconstruction
#' 
#' Wrapper function that returns just the optimal delay tau, suitable
#' for direct use in the Rosenstein Lyapunov exponent calculation.
#'
#' @param signal Numeric vector (time series)
#' @param n_bins Number of bins (default 64, as used in ACIER)
#' @param max_lag Maximum lag to search (default 100)
#' @return Integer: optimal delay tau (first minimum of AMI)
#' @examples
#' tau <- estimate_delay(acceleration_signal)
#' # Use tau in Rosenstein algorithm
#' 
#' @export
estimate_delay <- function(signal, n_bins = 64, max_lag = 100) {
  result <- ami(signal, n_bins = n_bins, n_lags = max_lag, plot = FALSE)
  
  if (is.na(result$optimal_lag)) {
    warning("No clear AMI minimum found. Using default delay of 10.")
    return(10L)
  }
  
  result$optimal_lag
}

