# ==============================================================================
# run_gait_analysis_batch.R
# ==============================================================================
# RACING Project - ORD 2025
# HE-Arc Health faculty, University of Applied Sciences Western Switzerland
#
# Master batch processing script for gait quality assessment from triaxial
# lumbar accelerometer data. Implements the complete ACIER study pipeline:
# tilt correction, signal standardization, and extraction of linear and
# nonlinear gait variability metrics.
#
# ## SCIENTIFIC BACKGROUND
#
# This pipeline computes validated gait quality metrics for fall risk assessment
# and gait complexity analysis. The methodology is based on:
#
# - **ACIER Study** (Attractor Complexity Index Empirical Rationalization):
#   Piergiovanni, S., & Terrier, P. (2024).
#   Effects of metronome walking on long-term attractor divergence and correlation
#   structure of gait: A validation study in older people.
#   Scientific Reports, 14, 15784.
#   https://doi.org/10.1038/s41598-024-65662-5
#
#   Piergiovanni, S., & Terrier, P. (2024).
#   Validity of linear and nonlinear measures of gait variability to characterize
#   aging gait with a single lower back accelerometer.
#   Sensors, 24(23), 7427.
#   https://doi.org/10.3390/s24237427
#
# - **Tilt Correction**:
#   Moe-Nilssen, R. (1998). A new method for evaluating motor control in gait
#   under real-life environmental conditions.
#   Clinical Biomechanics, 13(4-5), 320-327.
#
# - **Rosenstein Algorithm** (Lyapunov exponents):
#   Rosenstein, M. T., Collins, J. J., & De Luca, C. J. (1993).
#   A practical method for calculating largest Lyapunov exponents from small data sets.
#   Physica D, 65(1-2), 117-134.
#
# ## COMPUTED METRICS
#
# **Linear Metrics** (from truncated signals at original sampling rate):
#   - Step Frequency (SF): Dominant walking frequency via FFT
#   - RMS Norm: Movement intensity (proxy for walking speed)
#   - RMS Ratio: Mediolateral stability index (ML/norm)
#   - Step Regularity: First ACF peak (Fisher z-transformed)
#   - Stride Regularity: Second ACF peak (Fisher z-transformed)
#
# **Nonlinear Metrics** (from resampled signals with standardized length):
#   - AMI (x, y, z, norm): Optimal embedding delay (samples)
#   - LDS (x, y, z, norm): Local Dynamic Stability (0-0.5 strides)
#   - ACI (x, y, z, norm): Attractor Complexity Index (5-12 strides)
#
# ## PIPELINE ARCHITECTURE
#
# Per-file processing order (matches ACIER MATLAB implementation exactly):
#   1. CSV Import           → Raw triaxial acceleration data
#   2. Tilt Correction      → Remove gravity effects (Moe-Nilssen 1998)
#   3. Step Frequency       → FFT on vertical axis
#   4. Truncation           → Extract exactly N steps at original rate
#   5. Resampling           → Standardize to 75 samples/step
#   6. RMS (truncated)      → Movement intensity + ML ratio
#   7. ACF (truncated norm) → Step and stride regularity
#   8. AMI (resampled)      → Embedding delay for each axis (parallel)
#   9. Divergence (resampled) → Rosenstein algorithm (parallel)
#  10. LDS/ACI Fitting      → Short-term and long-term exponents (parallel)
#
# ## USAGE
#
# **Option 1: Command line**
#   Rscript run_gait_analysis_batch.R path/to/your_config.xlsx
#
# **Option 2: RStudio**
#   config_path <- "path/to/your_config.xlsx"
#   source("run_gait_analysis_batch.R")
#
# **Option 3: With custom settings**
#   config_path <- "my_config.xlsx"
#   verbose_batch <- FALSE        # Suppress per-step output
#   use_parallel <- FALSE         # Force sequential (for debugging)
#   scripts_folder <- "./scripts" # Custom script location
#   source("run_gait_analysis_batch.R")
#
# ## CONFIGURATION FILE
#
# The Excel configuration file (RACING_config_template.xlsx) must specify:
#   - data_folder: Path to CSV files
#   - output_folder: Where to write results
#   - study_name: Prefix for output filename
#   - All analysis parameters (see read_config.R documentation)
#
# ## OUTPUT
#
# **Primary output**: {study_name}_results.xlsx
#   - Sheet 1: Computed metrics for all files
#   - Sheet 2: Configuration used
#   - Sheet 3: Processing summary (success/warnings/failures)
#
# **Optional exports** (if enabled in config):
#   - divergence_curves/{filename}_{axis}_div.csv
#   - resampled_signals/{filename}_resampled.csv
#
# ## FEATURES
#
# - **Fault-tolerant**: One failed file doesn't stop the batch
# - **Parallel processing**: AMI + divergence computed across 4 axes simultaneously
# - **Progress tracking**: Per-step timing with ETA estimation
# - **Quality assurance**: Comprehensive warning collection and reporting
# - **Dependency checking**: Validates all required packages before starting
# - **Graceful degradation**: Optional packages (RANN, data.table) auto-detected
#
# ## DEPENDENCIES
#
# **Required packages**:
#   - openxlsx    (Excel I/O)
#   - readxl      (Config reading)
#   - gsignal OR signal (Resampling - gsignal preferred)
#
# **Recommended packages** (auto-detected, significant speedup):
#   - RANN        (Fast nearest neighbor search, ~10x speedup)
#   - data.table  (Fast CSV reading)
#   - parallel    (Multi-core processing - base R, usually available)
#
# Install with:
#   install.packages(c("openxlsx", "readxl", "gsignal", "RANN", "data.table"))
#
# ## PERFORMANCE
#
# **Typical processing times** (Intel i7, 4 cores):
#   - Short signals (250 steps): ~10-15 seconds/file
#   - Long signals (500 steps): ~20-30 seconds/file
#   - Speedup with parallel: ~2.5-3x vs sequential
#   - With RANN: ~1.5-2x additional speedup
#
# ## AUTHOR & LICENSE
#
# Author:  Philippe Terrier, PhD
#          HE-Arc Health, HES-SO
#          philippe.terrier@he-arc.ch
#
# Project: RACING
#          ORD 2025 Grant, HES-SO
#
# Date:    February 2026
#
# License: To be determined upon publication
#          (Likely GPL-3 or MIT for maximum scientific reusability)
#
# Citation: 
#
# ==============================================================================

# ==============================================================================
# 0. SETUP: Paths, Options, and Script Sourcing
# ==============================================================================

# -- 0a. Determine config file path -------------------------------------------
# Priority: (1) command-line argument, (2) variable set before source()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  config_path <- args[1]
} else if (!exists("config_path")) {
  # Default: look in working directory
  config_path <- "RACING_config_template.xlsx"
}

if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path,
       "\nUsage: Rscript run_gait_analysis_batch.R <config.xlsx>")
}

# -- 0b. User-configurable options --------------------------------------------
# Set these BEFORE source("run_gait_analysis_batch.R") to override defaults.

# Verbose: show per-step progress within each file
if (!exists("verbose_batch")) {
  verbose_batch <- TRUE
}

# Parallel: use multiple cores for AMI + divergence (4 axes)
if (!exists("use_parallel")) {
  use_parallel <- TRUE   # Set FALSE to force sequential
}

# -- 0c. Determine script directory --------------------------------------------
# All helper R scripts must be in the same folder as this batch script.
script_dir <- tryCatch({
  script_path <- sys.frame(1)$ofile
  if (!is.null(script_path)) dirname(normalizePath(script_path))
  else getwd()
}, error = function(e) getwd())

# Allow override
if (exists("scripts_folder")) {
  script_dir <- scripts_folder
}

# -- 0d. Source all required analysis scripts ----------------------------------
required_scripts <- c(
  "read_config.R",
  "import_csv_lowback.R",
  "tilt_correction.R",
  "step_frequency_detection.R",
  "resample_signal.R",
  "truncate_resample.R",
  "compute_rms.R",
  "acf_gait_regularity.R",
  "ami.R",
  "rosenstein_divergence.R",
  "fit_divergence_exponents.R",
  "write_results.R"
)


cat("\n")
cat("================================================================\n")
cat("  RACING Gait Analysis - Batch Processing\n")
cat("================================================================\n\n")

for (script_name in required_scripts) {
  script_file <- file.path(script_dir, script_name)
  if (!file.exists(script_file)) {
    stop(sprintf("Required script not found: %s\n  Expected in: %s",
                 script_name, script_dir))
  }
  source(script_file, local = FALSE)
}
cat(sprintf("  Loaded %d analysis scripts from: %s\n",
            length(required_scripts), script_dir))


# ==============================================================================
# 1. READ AND VALIDATE CONFIGURATION
# ==============================================================================

cat(sprintf("  Config file: %s\n\n", config_path))
config <- read_racing_config(config_path, verbose = TRUE)


# ==============================================================================
# 2. CHECK DEPENDENCIES
# ==============================================================================

pkg_check <- function(pkg, required = TRUE) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok && required) {
    stop(sprintf("Required package '%s' not installed. Run: install.packages('%s')",
                 pkg, pkg))
  }
  return(ok)
}

pkg_check("openxlsx", required = TRUE)
pkg_check("readxl", required = TRUE)

has_gsignal <- pkg_check("gsignal", required = FALSE)
has_signal  <- pkg_check("signal", required = FALSE)
if (!has_gsignal && !has_signal) {
  stop("Either 'gsignal' or 'signal' package is required for resampling.\n",
       "  Install with: install.packages('gsignal')  # preferred\n",
       "  Or fallback:  install.packages('signal')")
}

has_rann <- pkg_check("RANN", required = FALSE)
has_dt   <- pkg_check("data.table", required = FALSE)

cat("  Package status:\n")
cat("    openxlsx:   OK\n")
cat("    readxl:     OK\n")
cat(sprintf("    gsignal:    %s\n", ifelse(has_gsignal, "OK (preferred resampling)", "not found")))
cat(sprintf("    signal:     %s\n", ifelse(has_signal, "OK (fallback resampling)", "not found")))
cat(sprintf("    RANN:       %s\n", ifelse(has_rann, "OK (fast NN search)", "not found - using brute force")))
cat(sprintf("    data.table: %s\n", ifelse(has_dt, "OK (fast CSV reading)", "not found - using base R")))


# ==============================================================================
# 2b. PARALLEL SETUP
# ==============================================================================

# Determine available physical cores (cap at 4: one per axis)
n_cores <- min(4L, parallel::detectCores(logical = FALSE))
if (is.na(n_cores)) n_cores <- 1L
use_parallel <- use_parallel && (n_cores >= 2L)

if (use_parallel) {
  cl <- parallel::makeCluster(n_cores)

  # Source the required scripts directly on each worker.
  # This is more robust than clusterExport, which requires listing
  # every helper function by name (fragile if scripts change).
  worker_scripts <- file.path(script_dir, c(
    "ami.R",
    "rosenstein_divergence.R",
    "fit_divergence_exponents.R"
  ))
  parallel::clusterExport(cl, "worker_scripts")
  parallel::clusterEvalQ(cl, {
    for (s in worker_scripts) source(s, local = FALSE)
  })

  # Load RANN on workers if available (for fast nearest neighbor search)
  if (has_rann) {
    parallel::clusterEvalQ(cl, library(RANN))
  }

  cat(sprintf("  Parallel:     ON (%d cores for AMI + divergence)\n", n_cores))
} else {
  cl <- NULL
  cat("  Parallel:     OFF (single core)\n")
}


# ==============================================================================
# 3. DISCOVER CSV FILES
# ==============================================================================

csv_files <- list.files(config$data_folder, pattern = "\\.csv$",
                         full.names = TRUE, ignore.case = TRUE)

if (length(csv_files) == 0) {
  stop(sprintf("No CSV files found in: %s", config$data_folder))
}

csv_files <- sort(csv_files)

cat(sprintf("\n  Found %d CSV files in: %s\n\n",
            length(csv_files), config$data_folder))


# ==============================================================================
# 4. INITIALIZE RESULTS TABLE
# ==============================================================================

col_names <- c("file_id", "subject_id", "n_steps", "SF_hz",
               "RMS_norm", "RMS_ratio_ml_norm", "step_reg", "stride_reg",
               "AMI_x", "AMI_y", "AMI_z", "AMI_norm",
               "LDS_x", "LDS_y", "LDS_z", "LDS_norm",
               "ACI_x", "ACI_y", "ACI_z", "ACI_norm",
               "warnings")

results <- data.frame(matrix(NA, nrow = length(csv_files), ncol = length(col_names)),
                       stringsAsFactors = FALSE)
colnames(results) <- col_names

results$file_id    <- NA_integer_
results$subject_id <- NA_character_
results$warnings   <- NA_character_


# ==============================================================================
# 5. PRECOMPUTE DERIVED PARAMETERS
# ==============================================================================

# ACF max lag: config provides samples, function expects seconds
# ACF max lag: already precomputed by read_racing_config()
acf_max_lag_sec <- config$acf_max_lag_sec

# LDS range: already computed in read_config.R
lds_steps <- config$lds_steps

# ACI range: already in strides, bundled as vector in read_config.R
aci_stride_range <- config$aci_stride_range

# Axes to analyze
axes <- c("x", "y", "z", "norm")

# Timing
batch_start <- Sys.time()


# ==============================================================================
# 6. DEFINE PER-AXIS COMPUTATION FUNCTION (for parallel use)
# ==============================================================================
# Computes AMI + Rosenstein divergence + exponent fitting for one axis.
# Called 4 times (x, y, z, norm), either sequentially or via parLapply.

compute_axis_divergence <- function(args) {
  ax                <- args$ax
  signal            <- args$signal
  ami_n_bins        <- args$ami_n_bins
  ami_n_lags        <- args$ami_n_lags
  embedding_dim     <- args$embedding_dim
  mean_period       <- args$mean_period
  divergence_length <- args$divergence_length
  samples_per_step  <- args$samples_per_step
  lds_steps         <- args$lds_steps
  aci_stride_range  <- args$aci_stride_range

  warnings_ax <- character(0)

  # --- F: AMI ---
  ami_result <- tryCatch(
    ami(signal, n_bins = ami_n_bins, n_lags = ami_n_lags),
    error = function(e) {
      warnings_ax <<- c(warnings_ax, sprintf("AMI_%s: %s", ax, e$message))
      list(optimal_lag = NA_integer_)
    }
  )

  tau <- ami_result$optimal_lag
  if (is.na(tau)) {
    tau <- 10L
    warnings_ax <- c(warnings_ax, sprintf("AMI_%s: no minimum, using tau=10", ax))
  }

  # --- G1: Rosenstein divergence curve ---
  div_result <- tryCatch(
    rosenstein_divergence(
      x           = signal,
      m           = embedding_dim,
      tau         = tau,
      mean_period = mean_period,
      max_iter    = divergence_length
    ),
    error = function(e) {
      warnings_ax <<- c(warnings_ax, sprintf("Divergence_%s: %s", ax, e$message))
      NULL
    }
  )

  if (is.null(div_result)) {
    return(list(ax = ax, tau = tau, LDS = NA_real_, ACI = NA_real_,
                divergence = NULL, warnings = warnings_ax))
  }

  # --- G2: Fit LDS and ACI exponents ---
  exponents <- tryCatch(
    fit_divergence_exponents(
      divergence       = div_result$divergence,
      samples_per_step = samples_per_step,
      lds_steps        = lds_steps,
      aci_stride_range = aci_stride_range
    ),
    error = function(e) {
      warnings_ax <<- c(warnings_ax, sprintf("Fit_%s: %s", ax, e$message))
      list(LDS = NA_real_, ACI = NA_real_)
    }
  )

  list(ax = ax, tau = tau,
       LDS = exponents$LDS, ACI = exponents$ACI,
       divergence = div_result$divergence,
       warnings = warnings_ax)
}


# ==============================================================================
# 7. MAIN PROCESSING LOOP
# ==============================================================================

cat("----------------------------------------------------------------\n")
cat("  Processing files...\n")
cat("----------------------------------------------------------------\n\n")

n_success  <- 0L
n_failed   <- 0L
n_warnings <- 0L

# Track timing for ETA estimation
file_times <- numeric(length(csv_files))

for (i in seq_along(csv_files)) {

  file_path <- csv_files[i]
  file_name <- tools::file_path_sans_ext(basename(file_path))
  file_warnings <- character(0)
  file_start <- Sys.time()

  # Always fill file_id and subject_id
  results$file_id[i]    <- i
  results$subject_id[i] <- file_name

  # -- File header --
  if (verbose_batch) {
    cat(sprintf("  [%3d/%d] %s\n", i, length(csv_files), file_name))
  }

  tryCatch({

    # -- A. Import CSV ---------------------------------------------------------
    if (verbose_batch) cat("         A. Import CSV...")
    step_start <- Sys.time()

    raw_data <- import_csv_lowback(file_path, config, verbose = FALSE)

    import_warns <- attr(raw_data, "import_warnings")
    if (!is.null(import_warns)) {
      file_warnings <- c(file_warnings, import_warns)
    }

    acc_x <- raw_data$ap
    acc_y <- raw_data$v
    acc_z <- raw_data$ml
    
    if (verbose_batch) {
      cat(sprintf(" %d samples (%.0fs recording) [%.1fs]\n",
                  nrow(raw_data),
                  nrow(raw_data) / config$sampling_freq,
                  as.numeric(difftime(Sys.time(), step_start, units = "secs"))))
    }

    # -- B. Step Frequency Detection -------------------------------------------
    if (verbose_batch) cat("         B. Step frequency...")
    step_start <- Sys.time()

    sf_result <- step_frequency_fft(
      acc_y,
      sampling_freq = config$sampling_freq,
      freq_range    = c(config$sf_freq_min, config$sf_freq_max)
    )

    if (!is.null(sf_result$quality) && !is.null(sf_result$quality$warnings)) {
      for (w in sf_result$quality$warnings) {
        file_warnings <- c(file_warnings, paste0("SF: ", w))
      }
    }

    if (verbose_batch) {
      cat(sprintf(" %.3f Hz (%.0f steps/min) [%.1fs]\n",
                  sf_result$step_freq_hz,
                  sf_result$step_freq_hz * 60,
                  as.numeric(difftime(Sys.time(), step_start, units = "secs"))))
    }

    # -- C. Truncate & Resample ------------------------------------------------
    if (verbose_batch) cat("         C. Truncate & resample...")
    step_start <- Sys.time()

    prep <- truncate_resample(
      acc_x, acc_y, acc_z,
      sampling_freq         = config$sampling_freq,
      target_steps          = config$target_steps,
      samples_per_step      = config$samples_per_step,
      apply_tilt_correction = config$apply_tilt,
      auto_orient           = config$auto_orient,
      step_freq_hz          = sf_result$step_freq_hz
    )

    # Collect preprocessing warnings
    if (length(prep$quality$warnings) > 0) {
      file_warnings <- c(file_warnings, prep$quality$warnings)
    }

    if (prep$metadata$signal_too_short) {
      file_warnings <- c(file_warnings,
        sprintf("Signal too short: %d steps available (target: %d)",
                round(prep$metadata$actual_steps), config$target_steps))
    }

    if (!is.null(prep$tilt_correction)) {
      if (prep$tilt_correction$orientation_corrected) {
        file_warnings <- c(file_warnings, "Tilt: inverted orientation corrected")
      }
      max_tilt <- max(abs(prep$tilt_correction$theta_ap_deg),
                      abs(prep$tilt_correction$theta_ml_deg))
      if (max_tilt > 30) {
        file_warnings <- c(file_warnings,
          sprintf("Tilt: extreme angle %.1f deg", max_tilt))
      }
    }

    if (verbose_batch) {
      cat(sprintf(" %d -> %d samples (%d steps) [%.1fs]\n",
                  prep$metadata$truncated_length,
                  prep$metadata$resampled_length,
                  round(prep$metadata$actual_steps),
                  as.numeric(difftime(Sys.time(), step_start, units = "secs"))))
    }

    # -- D. RMS (on truncated signals) -----------------------------------------
    if (verbose_batch) cat("         D. RMS...")
    step_start <- Sys.time()

    rms <- compute_rms(
      az          = prep$truncated$z,
      norm_signal = prep$truncated$norm
    )

    if (verbose_batch) {
      cat(sprintf(" norm=%.4f ratio=%.3f [%.1fs]\n",
                  rms$rms_norm, rms$rms_ratio,
                  as.numeric(difftime(Sys.time(), step_start, units = "secs"))))
    }

    # -- E. ACF Gait Regularity (on truncated norm) ----------------------------
    if (verbose_batch) cat("         E. ACF regularity...")
    step_start <- Sys.time()

    reg <- tryCatch({
      compute_gait_regularity(
        acc_norm         = prep$truncated$norm,
        sampling_freq    = config$sampling_freq,
        max_lag_seconds  = acf_max_lag_sec,
        peak_order       = config$peak_min_distance,
        apply_fisher     = TRUE
      )
    }, error = function(e) {
      file_warnings <<- c(file_warnings, paste0("ACF: ", e$message))
      list(step_regularity = NA_real_, stride_regularity = NA_real_)
    })

    if (verbose_batch) {
      step_str <- if (is.na(reg$step_regularity)) "NA" else sprintf("%.4f", reg$step_regularity)
      stride_str <- if (is.na(reg$stride_regularity)) "NA" else sprintf("%.4f", reg$stride_regularity)
      cat(sprintf(" step=%s stride=%s [%.1fs]\n",
                  step_str, stride_str,
                  as.numeric(difftime(Sys.time(), step_start, units = "secs"))))
    }

    # -- F+G. AMI + Divergence + Fitting (4 axes, parallel) --------------------
    if (verbose_batch) {
      cat(sprintf("         F+G. AMI + Divergence (%s)...",
                  ifelse(use_parallel, sprintf("%d cores", n_cores), "sequential")))
    }
    step_start <- Sys.time()

    # Build argument list for each axis
    axis_args <- lapply(axes, function(ax) {
      list(ax                = ax,
           signal            = prep$resampled[[ax]],
           ami_n_bins        = config$ami_n_bins,
           ami_n_lags        = config$ami_n_lags,
           embedding_dim     = config$embedding_dim,
           mean_period       = config$mean_period,
           divergence_length = config$divergence_length,
           samples_per_step  = config$samples_per_step,
           lds_steps         = lds_steps,
           aci_stride_range  = aci_stride_range)
    })

    # Execute: parallel or sequential
    if (use_parallel) {
      axis_results <- parallel::parLapply(cl, axis_args, compute_axis_divergence)
    } else {
      axis_results <- lapply(axis_args, compute_axis_divergence)
    }

    fg_elapsed <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
    if (verbose_batch) cat(sprintf(" %.1fs\n", fg_elapsed))

    # Unpack results into vectors
    ami_vals <- setNames(integer(4), axes)
    lds_vals <- setNames(numeric(4), axes)
    aci_vals <- setNames(numeric(4), axes)
    div_curves <- list()  # Store for optional export (avoids recomputing)

    for (res in axis_results) {
      ax <- res$ax
      ami_vals[ax] <- res$tau
      lds_vals[ax] <- res$LDS
      aci_vals[ax] <- res$ACI
      file_warnings <- c(file_warnings, res$warnings)
      if (!is.null(res$divergence)) {
        div_curves[[ax]] <- res$divergence
      }
    }

    # Print per-axis summary
    if (verbose_batch) {
      for (res in axis_results) {
        lds_str <- if (is.na(res$LDS)) "   NA  " else sprintf("%.4f", res$LDS)
        aci_str <- if (is.na(res$ACI)) "   NA  " else sprintf("%.4f", res$ACI)
        cat(sprintf("              %4s: tau=%2d  LDS=%s  ACI=%s\n",
                    res$ax, res$tau, lds_str, aci_str))
      }
    }

    # -- H. Optional: Export divergence curves ---------------------------------
    # Uses curves already computed in G (no recomputation needed)
    if (isTRUE(config$export_div_curves)) {
      div_dir <- file.path(config$output_folder, "divergence_curves")
      if (!dir.exists(div_dir)) dir.create(div_dir, recursive = TRUE)

      for (ax in axes) {
        if (!is.null(div_curves[[ax]])) {
          div_df <- data.frame(
            sample = seq_along(div_curves[[ax]]),
            divergence = div_curves[[ax]]
          )
          write.csv(div_df,
                     file.path(div_dir, sprintf("%s_%s_div.csv", file_name, ax)),
                     row.names = FALSE)
        }
      }
    }

    # -- I. Optional: Export resampled signals ----------------------------------
    if (isTRUE(config$export_resampled)) {
      resamp_dir <- file.path(config$output_folder, "resampled_signals")
      if (!dir.exists(resamp_dir)) dir.create(resamp_dir, recursive = TRUE)

      resamp_df <- data.frame(
        x    = prep$resampled$x,
        y    = prep$resampled$y,
        z    = prep$resampled$z,
        norm = prep$resampled$norm
      )
      write.csv(resamp_df,
                 file.path(resamp_dir, sprintf("%s_resampled.csv", file_name)),
                 row.names = FALSE)
    }

    # -- J. Compile results row ------------------------------------------------
    results$n_steps[i]          <- prep$metadata$actual_steps
    results$SF_hz[i]            <- sf_result$step_freq_hz
    results$RMS_norm[i]         <- rms$rms_norm
    results$RMS_ratio_ml_norm[i] <- rms$rms_ratio
    results$step_reg[i]         <- reg$step_regularity
    results$stride_reg[i]       <- reg$stride_regularity
    results$AMI_x[i]            <- ami_vals["x"]
    results$AMI_y[i]            <- ami_vals["y"]
    results$AMI_z[i]            <- ami_vals["z"]
    results$AMI_norm[i]         <- ami_vals["norm"]
    results$LDS_x[i]            <- lds_vals["x"]
    results$LDS_y[i]            <- lds_vals["y"]
    results$LDS_z[i]            <- lds_vals["z"]
    results$LDS_norm[i]         <- lds_vals["norm"]
    results$ACI_x[i]            <- aci_vals["x"]
    results$ACI_y[i]            <- aci_vals["y"]
    results$ACI_z[i]            <- aci_vals["z"]
    results$ACI_norm[i]         <- aci_vals["norm"]

    # Join warnings
    results$warnings[i] <- if (length(file_warnings) > 0) {
      paste(file_warnings, collapse = "; ")
    } else {
      ""
    }

    n_success <- n_success + 1L
    if (length(file_warnings) > 0) n_warnings <- n_warnings + 1L

    # -- File summary with ETA -------------------------------------------------
    file_elapsed <- as.numeric(difftime(Sys.time(), file_start, units = "secs"))
    file_times[i] <- file_elapsed
    warn_flag <- if (length(file_warnings) > 0) {
      sprintf(" [%d warning(s)]", length(file_warnings))
    } else ""

    # Estimate time remaining
    avg_time <- mean(file_times[1:i])
    remaining <- (length(csv_files) - i) * avg_time
    eta_str <- if (remaining > 3600) {
      sprintf("%.1f h", remaining / 3600)
    } else if (remaining > 60) {
      sprintf("%.0f min", remaining / 60)
    } else {
      sprintf("%.0f s", remaining)
    }

    cat(sprintf("         -> DONE (%.1fs)%s  | ETA: %s remaining\n\n",
                file_elapsed, warn_flag, eta_str))

  }, error = function(e) {

    # -- Error handler: fill row with NAs --------------------------------------
    results$warnings[i] <<- paste0("ERROR: ", conditionMessage(e))
    n_failed <<- n_failed + 1L

    file_elapsed <- as.numeric(difftime(Sys.time(), file_start, units = "secs"))
    file_times[i] <<- file_elapsed

    cat(sprintf("         -> FAILED (%.1fs): %s\n\n",
                file_elapsed, conditionMessage(e)))
  })
}


# ==============================================================================
# 8. CLEAN UP PARALLEL CLUSTER
# ==============================================================================

if (use_parallel && !is.null(cl)) {
  parallel::stopCluster(cl)
  cat("  Parallel cluster stopped.\n")
}


# ==============================================================================
# 9. WRITE OUTPUT EXCEL FILE
# ==============================================================================

cat("\n----------------------------------------------------------------\n")
cat("  Writing results...\n")
cat("----------------------------------------------------------------\n\n")

if (!dir.exists(config$output_folder)) {
  dir.create(config$output_folder, recursive = TRUE)
}

output_file <- file.path(config$output_folder,
                          paste0(config$study_name, "_results.xlsx"))

write_racing_results(results, config, output_file)


# ==============================================================================
# 10. BATCH SUMMARY
# ==============================================================================

batch_elapsed <- as.numeric(difftime(Sys.time(), batch_start, units = "secs"))

cat("================================================================\n")
cat("  BATCH COMPLETE\n")
cat("================================================================\n")
cat(sprintf("  Files processed: %d / %d\n", n_success + n_failed, length(csv_files)))
cat(sprintf("  Successful:      %d\n", n_success))
cat(sprintf("  With warnings:   %d\n", n_warnings))
cat(sprintf("  Failed:          %d\n", n_failed))

# Format total time nicely
if (batch_elapsed > 3600) {
  cat(sprintf("  Total time:      %.1f min (%.1f h)\n",
              batch_elapsed / 60, batch_elapsed / 3600))
} else if (batch_elapsed > 60) {
  cat(sprintf("  Total time:      %.1f min\n", batch_elapsed / 60))
} else {
  cat(sprintf("  Total time:      %.1f seconds\n", batch_elapsed))
}
cat(sprintf("  Average:         %.1f s/file\n", batch_elapsed / length(csv_files)))
if (use_parallel) {
  cat(sprintf("  Parallelism:     %d cores used\n", n_cores))
}
cat(sprintf("  Output:          %s\n", output_file))
cat("================================================================\n\n")

# List failed files if any
if (n_failed > 0) {
  failed_idx <- which(grepl("^ERROR:", results$warnings))
  cat("  FAILED FILES:\n")
  for (idx in failed_idx) {
    cat(sprintf("    %s: %s\n", results$subject_id[idx], results$warnings[idx]))
  }
  cat("\n")
}

# List files with warnings if any
if (n_warnings > 0) {
  warn_idx <- which(!is.na(results$warnings) & nchar(results$warnings) > 0 &
                      !grepl("^ERROR:", results$warnings))
  if (length(warn_idx) > 0) {
    cat("  FILES WITH WARNINGS:\n")
    for (idx in warn_idx) {
      cat(sprintf("    %s: %s\n", results$subject_id[idx], results$warnings[idx]))
    }
    cat("\n")
  }
}
