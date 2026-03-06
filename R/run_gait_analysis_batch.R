# ==============================================================================
# run_gait_analysis_batch.R
# ==============================================================================
# RACING Project - ORD 2025
# Master batch script: Reads configuration, processes all accelerometer CSV
# files, computes gait quality metrics, and writes results to Excel.
#
# Usage:
#   1. Edit RACING_config_template.xlsx with your paths and parameters
#   2. Run this script:
#      Rscript run_gait_analysis_batch.R path/to/your_config.xlsx
#      -- OR --
#      In RStudio: set config_path below and source() this file
#
# Pipeline per file:
#   CSV  import  tilt correction  truncation  resampling 
#   RMS + ACF regularity (on truncated) 
#   AMI  Rosenstein divergence  LDS/ACI fitting (on resampled) 
#   results row
#
# Parallelization:
#   The nonlinear analysis (AMI + Rosenstein + LDS/ACI fitting) is the
#   computational bottleneck (~90% of per-file processing time). The four
#   axes (AP, V, ML, norm) are independent and processed in parallel across
#   cluster workers using the parallel package. A PSOCK cluster is created
#   at startup (default: 4 workers, one per axis) and reused for all files.
#   Falls back to sequential lapply when parallel setup fails or on
#   single-core systems.
#
# Safety features:
#   - All output paths resolved to absolute at startup
#   - Incremental CSV checkpoint every N files (default: 50)
#   - RDS backup saved before Excel generation
#   - Post-save file verification
#   - CSV fallback if Excel save fails
#
# Dependencies:
#   Required:    openxlsx, gsignal (or signal as fallback)
#   Recommended: RANN (10x speedup for nearest neighbor search)
#   Optional:    data.table (faster CSV reading)
#   Base R:      parallel (used for across-axis parallelization)
#
# Author: Philippe Terrier / HE-Arc
# Date:   February 2026
# ==============================================================================


# ==============================================================================
# 0. SETUP: Paths and Script Sourcing
# ==============================================================================

#  0a. Determine config file path 
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

#  0b. Determine script directory 
# All helper R scripts must be in the same folder as this batch script.
# Auto-detect: if running via Rscript, use the script's own location;
# otherwise, fall back to working directory.
script_dir <- tryCatch({
  # Works when called via Rscript or source()
  script_path <- sys.frame(1)$ofile
  if (!is.null(script_path)) dirname(normalizePath(script_path))
  else getwd()
}, error = function(e) getwd())

# Allow override: if user sets `scripts_folder` before sourcing, use that
if (exists("scripts_folder")) {
  script_dir <- scripts_folder
}

#  0c. Source all required analysis scripts 
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
# 1b. RESOLVE ALL OUTPUT PATHS TO ABSOLUTE (prevent ghost-file problem)
# ==============================================================================

config$output_folder <- normalizePath(config$output_folder, mustWork = FALSE)
cat(sprintf("\n  Output folder (absolute): %s\n", config$output_folder))
cat(sprintf("  Working directory:        %s\n\n", getwd()))


# ==============================================================================
# 2. CHECK DEPENDENCIES
# ==============================================================================

# Required
pkg_check <- function(pkg, required = TRUE) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok && required) {
    stop(sprintf("Required package '%s' not installed. Run: install.packages('%s')",
                 pkg, pkg))
  }
  return(ok)
}

pkg_check("openxlsx", required = TRUE)

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
cat(sprintf("    openxlsx:   OK\n"))
cat(sprintf("    gsignal:    %s\n", ifelse(has_gsignal, "OK (preferred resampling)", "not found")))
cat(sprintf("    signal:     %s\n", ifelse(has_signal, "OK (fallback resampling)", "not found")))
cat(sprintf("    RANN:       %s\n", ifelse(has_rann, "OK (fast NN search)", "not found - using brute force")))
cat(sprintf("    data.table: %s\n", ifelse(has_dt, "OK (fast CSV reading)", "not found - using base R")))


# ==============================================================================
# 2b. PARALLEL CLUSTER SETUP (across-axis parallelization)
# ==============================================================================
# The four axes (x, y, z, norm) are independent for AMI, Rosenstein divergence,
# and LDS/ACI fitting. We create a PSOCK cluster with up to 4 workers (one per
# axis) and source the three required scripts on each worker. The cluster is
# reused for all files in the batch. If setup fails, processing falls back to
# sequential lapply with identical results.

N_PARALLEL_WORKERS <- 4L   # one per axis (x, y, z, norm)

cl <- NULL
use_parallel <- FALSE

tryCatch({
  n_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(n_cores)) n_cores <- 1L
  n_workers <- min(N_PARALLEL_WORKERS, max(1L, n_cores))

  if (n_workers > 1L) {
    cl <- parallel::makeCluster(n_workers)

    # Source analysis scripts on each worker (project convention:
    # use clusterEvalQ(source(...)) rather than clusterExport for functions)
    parallel::clusterExport(cl, "script_dir", envir = environment())
    parallel::clusterEvalQ(cl, {
      source(file.path(script_dir, "ami.R"))
      source(file.path(script_dir, "rosenstein_divergence.R"))
      source(file.path(script_dir, "fit_divergence_exponents.R"))
    })

    use_parallel <- TRUE
    cat(sprintf("\n  Parallel:     %d workers (PSOCK cluster, %d cores detected)\n",
                n_workers, n_cores))
  } else {
    cat("\n  Parallel:     disabled (single core detected)\n")
  }
}, error = function(e) {
  cat(sprintf("\n  Parallel:     disabled (%s)\n", e$message))
  cl <<- NULL
  use_parallel <<- FALSE
})


# ==============================================================================
# 3. DISCOVER CSV FILES
# ==============================================================================

csv_files <- list.files(config$data_folder, pattern = "\\.csv$",
                         full.names = TRUE, ignore.case = TRUE)

if (length(csv_files) == 0) {
  stop(sprintf("No CSV files found in: %s", config$data_folder))
}

# Sort for reproducibility
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

# Pre-allocate with NA
results <- data.frame(matrix(NA, nrow = length(csv_files), ncol = length(col_names)),
                       stringsAsFactors = FALSE)
colnames(results) <- col_names

# Force correct types
results$file_id    <- NA_integer_
results$subject_id <- NA_character_
results$warnings   <- NA_character_


# ==============================================================================
# 5. PRECOMPUTE DERIVED PARAMETERS
# ==============================================================================

# ACF max lag: config provides samples, function expects seconds
acf_max_lag_sec <- config$acf_max_lag / config$sampling_freq

# LDS range: config provides strides (0.5), function expects steps (1)
lds_steps <- config$lds_steps  # already computed in read_config.R

# ACI range: already in strides, bundled as vector in read_config.R
aci_stride_range <- config$aci_stride_range

# Axes to analyze for AMI / divergence / fitting
axes <- c("x", "y", "z", "norm")

# Timing
batch_start <- Sys.time()

# -- Checkpoint settings --
CHECKPOINT_INTERVAL <- 50L   # Save CSV checkpoint every N files

# -- Export analysis parameters to parallel workers --
# These scalar values are captured in the worker function closure and sent
# to each worker via parLapply serialization. We store them with a .cfg_
# prefix in the batch script scope to keep them accessible to the closure.
.cfg_ami_bins   <- config$ami_n_bins
.cfg_ami_lags   <- config$ami_n_lags
.cfg_emb_dim    <- config$embedding_dim
.cfg_mean_period <- config$mean_period
.cfg_div_length <- config$divergence_length
.cfg_sps        <- config$samples_per_step
.cfg_lds_steps  <- lds_steps
.cfg_aci_range  <- aci_stride_range


# ==============================================================================
# 6. MAIN PROCESSING LOOP
# ==============================================================================

cat("\n")
cat("  Processing files...\n")
cat("\n\n")

# -- Ensure output folder exists before any file writing --
if (!dir.exists(config$output_folder)) {
  dir.create(config$output_folder, recursive = TRUE)
  cat(sprintf("  Created output folder: %s\n", config$output_folder))
}

# -- Verify output folder is writable --
.write_test <- file.path(config$output_folder, ".racing_write_test")
.write_ok <- tryCatch({
  writeLines("test", .write_test)
  file.remove(.write_test)
  TRUE
}, error = function(e) FALSE)

if (!.write_ok) {
  stop("Output folder is NOT writable: ", config$output_folder,
       "\n  Check permissions or choose a different path.")
}

# -- Debug resampling log (set to NULL to disable) --
resample_log_path <- file.path(config$output_folder, "resample_log.txt")
resample_log_header(resample_log_path)
cat(sprintf("  Resampling debug log: %s\n", resample_log_path))

# -- Checkpoint paths --
checkpoint_csv <- file.path(config$output_folder,
                             paste0(config$study_name, "_checkpoint.csv"))
checkpoint_rds <- file.path(config$output_folder,
                             paste0(config$study_name, "_checkpoint.rds"))
cat(sprintf("  Checkpoint every:    %d files -> %s\n\n", 
            CHECKPOINT_INTERVAL, checkpoint_csv))

n_success  <- 0L
n_failed   <- 0L
n_warnings <- 0L

for (i in seq_along(csv_files)) {

  file_path <- csv_files[i]
  file_name <- tools::file_path_sans_ext(basename(file_path))
  file_warnings <- character(0)
  file_start <- Sys.time()

  # Always fill file_id and subject_id
  results$file_id[i]    <- i
  results$subject_id[i] <- file_name

  tryCatch({

    #  A. Import CSV 
    raw_data <- import_csv_lowback(file_path, config, verbose = FALSE)

    # Collect any import warnings
    import_warns <- attr(raw_data, "import_warnings")
    if (!is.null(import_warns)) {
      file_warnings <- c(file_warnings, import_warns)
    }

    acc_x <- raw_data$ap
    acc_y <- raw_data$v
    acc_z <- raw_data$ml

    # -- B. Truncate & Resample (includes tilt correction + SF detection) --
    # SF is detected INTERNALLY by truncate_resample() on the tilt-corrected
    # vertical signal, matching the MATLAB pipeline order:
    #   raw -> tilt_correction -> SF detection -> truncate -> resample
    # The SF value used is returned in prep$metadata$step_freq_hz.
    
    resample_log_sep(resample_log_path,
                     sprintf("FILE [%d/%d]: %s", i, length(csv_files), file_name))
    
    prep <- truncate_resample(
      acc_x, acc_y, acc_z,
      sampling_freq        = config$sampling_freq,
      target_steps         = config$target_steps,
      samples_per_step     = config$samples_per_step,
      apply_tilt_correction = config$apply_tilt,
      auto_orient          = config$auto_orient,
      freq_range           = c(config$sf_freq_min, config$sf_freq_max),
      resampling_method    = "direct",
      log_file             = resample_log_path
    )

    # Collect preprocessing warnings
    if (length(prep$quality$warnings) > 0) {
      file_warnings <- c(file_warnings, prep$quality$warnings)
    }

    # Signal-too-short warning
    if (prep$metadata$signal_too_short) {
      file_warnings <- c(file_warnings,
        sprintf("Signal too short: %d steps available (target: %d)",
                round(prep$metadata$actual_steps), config$target_steps))
    }

    # Tilt correction diagnostics
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

    #  D. RMS (on truncated signals) 
    rms <- compute_rms(
      az          = prep$truncated$z,
      norm_signal = prep$truncated$norm
    )

    #  E. ACF Gait Regularity (on truncated norm) 
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

    # ---- F+G. Nonlinear metrics: AMI + Rosenstein + LDS/ACI (4 axes) --------
    # Each axis is independent: AMI determines the optimal embedding delay
    # (tau), Rosenstein's algorithm computes the divergence curve, and linear
    # regression extracts LDS and ACI. When a parallel cluster is available,
    # the four axes run simultaneously (one per worker); otherwise they are
    # processed sequentially with identical results.

    ami_vals <- setNames(integer(4), axes)
    lds_vals <- setNames(numeric(4), axes)
    aci_vals <- setNames(numeric(4), axes)

    # Capture per-file resampled signals for the worker closure
    resampled_signals <- prep$resampled

    # Worker function: full nonlinear pipeline for one axis.
    # References resampled_signals and .cfg_* variables from enclosing scope;
    # these are serialized into the closure by parLapply and accessible via
    # lexical scoping by lapply.
    .run_axis <- function(ax) {
      ax_warn <- character(0)

      # F: AMI -- optimal embedding delay
      ami_result <- tryCatch(
        ami(resampled_signals[[ax]],
            n_bins = .cfg_ami_bins,
            n_lags = .cfg_ami_lags),
        error = function(e) {
          ax_warn <<- c(ax_warn, sprintf("AMI_%s: %s", ax, e$message))
          list(optimal_lag = NA_integer_)
        }
      )

      tau <- ami_result$optimal_lag
      if (is.na(tau)) {
        tau <- 10L
        ax_warn <- c(ax_warn,
          sprintf("AMI_%s: no minimum found, using tau=10", ax))
      }

      # G1: Rosenstein divergence curve
      div_result <- tryCatch(
        rosenstein_divergence(
          x           = resampled_signals[[ax]],
          m           = .cfg_emb_dim,
          tau         = tau,
          mean_period = .cfg_mean_period,
          max_iter    = .cfg_div_length
        ),
        error = function(e) {
          ax_warn <<- c(ax_warn, sprintf("Divergence_%s: %s", ax, e$message))
          NULL
        }
      )

      lds_val <- NA_real_
      aci_val <- NA_real_

      if (!is.null(div_result)) {
        # G2: Fit LDS and ACI exponents
        exponents <- tryCatch(
          fit_divergence_exponents(
            divergence       = div_result$divergence,
            samples_per_step = .cfg_sps,
            lds_steps        = .cfg_lds_steps,
            aci_stride_range = .cfg_aci_range
          ),
          error = function(e) {
            ax_warn <<- c(ax_warn, sprintf("Fit_%s: %s", ax, e$message))
            list(LDS = NA_real_, ACI = NA_real_)
          }
        )
        lds_val <- exponents$LDS
        aci_val <- exponents$ACI
      }

      list(axis = ax, ami = tau, lds = lds_val, aci = aci_val,
           warnings = ax_warn)
    }

    # Dispatch: parallel (parLapply) or sequential (lapply)
    if (use_parallel) {
      axis_results <- parallel::parLapply(cl, axes, .run_axis)
    } else {
      axis_results <- lapply(axes, .run_axis)
    }

    # Collect results from all axes
    for (res in axis_results) {
      ami_vals[res$axis] <- res$ami
      lds_vals[res$axis] <- res$lds
      aci_vals[res$axis] <- res$aci
      file_warnings <- c(file_warnings, res$warnings)
    }

    #  H. Optional: Export divergence curves 
    if (isTRUE(config$export_div_curves)) {
      div_dir <- file.path(config$output_folder, "divergence_curves")
      if (!dir.exists(div_dir)) dir.create(div_dir, recursive = TRUE)

      for (ax in axes) {
        div_result <- tryCatch({
          rosenstein_divergence(
            x           = prep$resampled[[ax]],
            m           = config$embedding_dim,
            tau         = ami_vals[ax],
            mean_period = config$mean_period,
            max_iter    = config$divergence_length
          )
        }, error = function(e) NULL)

        if (!is.null(div_result)) {
          div_df <- data.frame(
            sample = seq_along(div_result$divergence),
            divergence = div_result$divergence
          )
          write.csv(div_df,
                     file.path(div_dir, sprintf("%s_%s_div.csv", file_name, ax)),
                     row.names = FALSE)
        }
      }
    }

    #  I. Optional: Export resampled signals 
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

    #  J. Compile results row 
    results$n_steps[i]          <- prep$metadata$actual_steps
    results$SF_hz[i]            <- prep$metadata$step_freq_hz
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

    elapsed <- round(as.numeric(difftime(Sys.time(), file_start, units = "secs")), 1)
    warn_flag <- if (length(file_warnings) > 0) " [!]" else ""
    cat(sprintf("  [%3d/%d] %-35s  OK  (%.1fs)%s\n",
                i, length(csv_files), file_name, elapsed, warn_flag))

  }, error = function(e) {

    #  Error handler: fill row with NAs 
    results$warnings[i] <<- paste0("ERROR: ", conditionMessage(e))
    n_failed <<- n_failed + 1L

    cat(sprintf("  [%3d/%d] %-35s  FAILED: %s\n",
                i, length(csv_files), file_name, conditionMessage(e)))
  })

  # -- Incremental checkpoint save --
  if (i %% CHECKPOINT_INTERVAL == 0 || i == length(csv_files)) {
    tryCatch({
      write.csv(results[1:i, ], checkpoint_csv, row.names = FALSE)
      saveRDS(results[1:i, ], checkpoint_rds)
    }, error = function(e) {
      cat(sprintf("  [!] Checkpoint save failed at file %d: %s\n", i, e$message))
    })
  }
}


# ==============================================================================
# 7. SHUTDOWN PARALLEL CLUSTER
# ==============================================================================

if (!is.null(cl)) {
  parallel::stopCluster(cl)
  cl <- NULL
}


# ==============================================================================
# 8. WRITE OUTPUT EXCEL FILE
# ==============================================================================

cat("\n\n")
cat("  Writing results...\n")
cat("\n\n")

# (output folder already ensured at start of processing loop)

output_file <- file.path(config$output_folder,
                          paste0(config$study_name, "_results.xlsx"))

write_racing_results(results, config, output_file)


# ==============================================================================
# 9. BATCH SUMMARY
# ==============================================================================

batch_elapsed <- as.numeric(difftime(Sys.time(), batch_start, units = "secs"))

# Resolve for display
output_file_abs <- normalizePath(output_file, mustWork = FALSE)

cat("================================================================\n")
cat("  BATCH COMPLETE\n")
cat("================================================================\n")
cat(sprintf("  Files processed: %d / %d\n", n_success, length(csv_files)))
cat(sprintf("  Successful:      %d\n", n_success))
cat(sprintf("  With warnings:   %d\n", n_warnings))
cat(sprintf("  Failed:          %d\n", n_failed))
cat(sprintf("  Total time:      %.1f seconds (%.1f s/file)\n",
            batch_elapsed,
            batch_elapsed / length(csv_files)))
cat(sprintf("  Parallelization: %s\n",
            ifelse(use_parallel,
                   sprintf("ON (%d workers, across axes)", N_PARALLEL_WORKERS),
                   "OFF (sequential)")))
cat(sprintf("  Output:          %s\n", output_file_abs))

# Verify output exists
if (file.exists(output_file_abs)) {
  cat(sprintf("  File size:       %.1f KB\n", file.size(output_file_abs) / 1024))
  cat("  Status:          [OK] File verified on disk\n")
} else {
  cat("  Status:          [!!] FILE NOT FOUND -- check checkpoint files\n")
  cat(sprintf("  Checkpoint CSV:  %s\n", checkpoint_csv))
  cat(sprintf("  Checkpoint RDS:  %s\n", checkpoint_rds))
}

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

# -- Clean up checkpoint files on success --
if (file.exists(output_file_abs) && file.size(output_file_abs) > 0) {
  # Keep RDS backup (small, useful), remove CSV checkpoint
  if (file.exists(checkpoint_csv)) {
    file.remove(checkpoint_csv)
  }
  cat("  Checkpoint CSV cleaned up. RDS backup retained.\n\n")
} else {
  cat("  [!] Checkpoint files RETAINED (output file missing or empty).\n")
  cat("  [!] To recover, run in R console:\n")
  cat(sprintf("      source('write_results.R')\n"))
  cat(sprintf("      results_df <- readRDS('%s')\n", checkpoint_rds))
  cat(sprintf("      write_racing_results(results_df, config, 'recovered.xlsx')\n\n"))
}
