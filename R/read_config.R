# read_config.R
# ==============
# RACING Project - ORD 2025
# Helper script: Read and validate the RACING configuration Excel file
#
# Usage:
#   source("read_config.R")
#   config <- read_racing_config("RACING_config.xlsx")
#
# Returns a named list with all configuration parameters, validated and
# type-converted, ready for use by the batch processing script.
#
# Dependencies: readxl (for reading config), openxlsx (for writing results)
# ==============

# Check that readxl is available (reads xlsx via C++ - no unzip issues on Windows)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install with: install.packages('readxl')")
}

#' Read RACING configuration from Excel file
#'
#' Reads all parameters from the '1_Configuration' sheet of a RACING config
#' Excel file, validates them, and returns a named list ready for batch use.
#'
#' @param xlsx_path Path to the RACING configuration Excel file
#' @param verbose   If TRUE, print a summary of loaded parameters (default TRUE)
#'
#' @return A named list with all configuration parameters, type-converted.
#'   Returns NULL (with error messages) if critical validation fails.
#'
#' @examples
#' config <- read_racing_config("RACING_config.xlsx")
#' config$sampling_freq   # 256
#' config$target_steps    # 250
#'
read_racing_config <- function(xlsx_path, verbose = TRUE) {
  
  # -- 1. Check file exists --------------------------------------------------
  if (!file.exists(xlsx_path)) {
    stop(sprintf("Configuration file not found: %s", xlsx_path))
  }
  
  # -- 2. Load the workbook --------------------------------------------------
  # Use readxl instead of openxlsx for reading: readxl uses a C++ library
  # (RapidXML) that reads xlsx directly without R's unzip(), which fails
  # on Windows with OneDrive/long/special-character paths.
  cfg_sheet <- as.data.frame(
    readxl::read_excel(xlsx_path,
                       sheet = "1_Configuration",
                       col_names = FALSE,
                       col_types = "text",
                       .name_repair = "minimal")
  )
  
  # -- 3. Helper to read a parameter from a specific Excel row ---------------
  # readxl::read_excel with col_names = FALSE preserves row positions
  # (row 1 of data frame = row 1 of Excel), so we read directly by row number.
  read_param <- function(excel_row) {
    val <- cfg_sheet[excel_row, 3]  # Column C = column 3
    if (is.null(val) || length(val) == 0) return(NA)
    return(val)
  }
  
  # Helper for type conversion
  to_integer <- function(val, param_name) {
    if (is.na(val)) return(NA_integer_)
    result <- suppressWarnings(as.integer(val))
    if (is.na(result)) {
      warning(sprintf("Could not convert '%s' to integer for parameter '%s'", val, param_name))
    }
    return(result)
  }
  
  to_numeric <- function(val, param_name) {
    if (is.na(val)) return(NA_real_)
    result <- suppressWarnings(as.numeric(val))
    if (is.na(result)) {
      warning(sprintf("Could not convert '%s' to numeric for parameter '%s'", val, param_name))
    }
    return(result)
  }
  
  to_logical <- function(val, param_name) {
    if (is.na(val)) return(NA)
    if (is.logical(val)) return(val)
    val_upper <- toupper(trimws(as.character(val)))
    if (val_upper %in% c("TRUE", "T", "YES", "1")) return(TRUE)
    if (val_upper %in% c("FALSE", "F", "NO", "0")) return(FALSE)
    warning(sprintf("Could not convert '%s' to logical for parameter '%s'", val, param_name))
    return(NA)
  }
  
  to_string <- function(val, param_name) {
    if (is.na(val)) return(NA_character_)
    return(trimws(as.character(val)))
  }
  
  # -- 4. Read all parameters by Excel row -----------------------------------
  # Row references from RACING_Batch_Script_Blueprint_v2.md, Section 1
  
  config <- list(
    # --- STEP 1: Input Data ---
    data_folder    = to_string(read_param(8),  "data_folder"),
    sampling_freq  = to_integer(read_param(9), "sampling_frequency"),
    col_ap         = to_integer(read_param(10), "column_ap"),
    col_v          = to_integer(read_param(11), "column_v"),
    col_ml         = to_integer(read_param(12), "column_ml"),
    accel_units    = to_string(read_param(13), "acceleration_units"),
    csv_header     = to_logical(read_param(14), "csv_has_header"),
    csv_delim      = to_string(read_param(15), "csv_delimiter"),
    
    # --- STEP 2: Preprocessing (Tilt Correction) ---
    apply_tilt     = to_logical(read_param(20), "apply_tilt_correction"),
    auto_orient    = to_logical(read_param(21), "auto_detect_orientation"),
    
    # --- STEP 3: Signal Truncation & Resampling ---
    target_steps     = to_integer(read_param(26), "target_steps"),
    samples_per_step = to_integer(read_param(27), "samples_per_step"),
    # C28 is a formula (=C26*C27); openxlsx reads the cached value if saved,
    # but we compute it explicitly for safety
    target_samples   = NA_integer_,  # computed below
    
    # --- STEP 4: Step Frequency Detection ---
    sf_freq_min    = to_numeric(read_param(33), "sf_freq_min"),
    sf_freq_max    = to_numeric(read_param(34), "sf_freq_max"),
    
    # --- STEP 5: Linear Metrics (RMS, Regularity) ---
    acf_max_lag       = to_integer(read_param(39), "acf_max_lag"),
    peak_min_distance = to_integer(read_param(40), "peak_min_distance"),
    
    # --- STEP 6: Phase Space Reconstruction (AMI) ---
    ami_n_bins     = to_integer(read_param(45), "ami_n_bins"),
    ami_n_lags     = to_integer(read_param(46), "ami_n_lags"),
    embedding_dim  = to_integer(read_param(47), "embedding_dimension"),
    
    # --- STEP 7: Divergence Analysis (LDS & ACI) ---
    divergence_length = to_integer(read_param(52), "divergence_length"),
    mean_period       = to_integer(read_param(53), "mean_period"),
    lds_range_end     = to_numeric(read_param(54), "lds_range_end"),
    aci_range_start   = to_numeric(read_param(55), "aci_range_start"),
    aci_range_end     = to_numeric(read_param(56), "aci_range_end"),
    
    # --- STEP 8: Output Options ---
    output_folder     = to_string(read_param(60), "output_folder"),
    study_name        = to_string(read_param(61), "study_name"),
    export_div_curves = to_logical(read_param(62), "export_divergence_curves"),
    export_resampled  = to_logical(read_param(63), "export_resampled_signals")
  )
  
  # -- 5. Compute derived parameters ----------------------------------------
  
  # target_samples: compute from components (don't rely on Excel formula cache)
  if (!is.na(config$target_steps) && !is.na(config$samples_per_step)) {
    config$target_samples <- config$target_steps * config$samples_per_step
  }
  
  # Precompute batch-useful derived values (see Blueprint Section 5)
  # ACF max lag in seconds (function expects seconds, config provides samples)
  if (!is.na(config$acf_max_lag) && !is.na(config$sampling_freq)) {
    config$acf_max_lag_sec <- config$acf_max_lag / config$sampling_freq
  } else {
    config$acf_max_lag_sec <- NA_real_
  }
  
  # LDS range in steps (function expects steps, config provides strides)
  # 0.5 strides x 2 = 1 step
  if (!is.na(config$lds_range_end)) {
    config$lds_steps <- as.integer(round(config$lds_range_end * 2))
  } else {
    config$lds_steps <- NA_integer_
  }
  
  # ACI range in strides (already in strides, just bundle as a vector)
  if (!is.na(config$aci_range_start) && !is.na(config$aci_range_end)) {
    config$aci_stride_range <- c(config$aci_range_start, config$aci_range_end)
  } else {
    config$aci_stride_range <- c(NA_real_, NA_real_)
  }
  
  # Normalize tab delimiter: Excel users may type "tab" or "TAB" as text,
  # but import_csv_lowback() needs the actual tab character "\t"
  if (!is.na(config$csv_delim) && toupper(trimws(config$csv_delim)) == "TAB") {
    config$csv_delim <- "\t"
  }
  
  # Store source file path for traceability
  config$config_file <- normalizePath(xlsx_path, mustWork = FALSE)
  
  # -- 6. Validate -----------------------------------------------------------
  validation <- validate_racing_config(config)
  
  if (length(validation$errors) > 0) {
    cat("\n+== CONFIGURATION ERRORS =========================================+\n")
    for (err in validation$errors) {
      cat(sprintf("  [X] %s\n", err))
    }
    cat("+================================================================+\n\n")
  }
  
  if (length(validation$warnings) > 0) {
    cat("\n+-- Configuration Warnings ---------------------------------------+\n")
    for (wrn in validation$warnings) {
      cat(sprintf("  [!] %s\n", wrn))
    }
    cat("+----------------------------------------------------------------+\n\n")
  }
  
  if (!validation$valid) {
    stop("Configuration validation failed. Fix the errors above before proceeding.")
  }
  
  # -- 7. Print summary -----------------------------------------------------
  if (verbose) {
    print_config_summary(config)
  }
  
  return(config)
}


#' Validate RACING configuration parameters
#'
#' Checks for missing required parameters, invalid types, and out-of-range values.
#'
#' @param config Named list of configuration parameters
#' @return List with: valid (logical), errors (character vector), warnings (character vector)
#'
validate_racing_config <- function(config) {
  
  errors   <- character(0)
  warnings <- character(0)
  
  # --- Required parameters (batch cannot run without these) ---
  required_params <- c("data_folder", "output_folder", "sampling_freq",
                       "col_ap", "col_v", "col_ml", "target_steps",
                       "samples_per_step", "study_name")
  
  for (param in required_params) {
    val <- config[[param]]
    if (is.null(val) || (length(val) == 1 && is.na(val))) {
      errors <- c(errors, sprintf("Missing required parameter: %s", param))
    }
  }
  
  # If critical params are missing, skip range checks (they'd fail on NA)
  if (length(errors) > 0) {
    return(list(valid = FALSE, errors = errors, warnings = warnings))
  }
  
  # --- Path validation ---
  if (!dir.exists(config$data_folder)) {
    errors <- c(errors, sprintf("Data folder not found: '%s'", config$data_folder))
  }
  
  # Output folder: create if it doesn't exist (warning, not error)
  if (!dir.exists(config$output_folder)) {
    warnings <- c(warnings,
                  sprintf("Output folder does not exist yet: '%s' (will be created at runtime)",
                          config$output_folder))
  }
  
  # --- Numeric range checks ---
  if (config$sampling_freq < 50 || config$sampling_freq > 1000) {
    warnings <- c(warnings,
                  sprintf("Unusual sampling frequency: %d Hz (typical: 100-256 Hz)",
                          config$sampling_freq))
  }
  
  if (config$target_steps < 50) {
    warnings <- c(warnings,
                  sprintf("target_steps = %d is very low; nonlinear metrics (LDS, ACI) need >= 200 steps for reliability",
                          config$target_steps))
  }
  
  if (config$samples_per_step < 20 || config$samples_per_step > 200) {
    warnings <- c(warnings,
                  sprintf("Unusual samples_per_step: %d (ACIER default: 75)",
                          config$samples_per_step))
  }
  
  # Column indices must be positive and distinct
  cols <- c(config$col_ap, config$col_v, config$col_ml)
  if (any(cols < 1)) {
    errors <- c(errors, "Column indices must be >= 1")
  }
  if (length(unique(cols)) != 3) {
    errors <- c(errors,
                sprintf("Column indices must be distinct (got AP=%d, V=%d, ML=%d)",
                        config$col_ap, config$col_v, config$col_ml))
  }
  
  # Acceleration units
  if (!config$accel_units %in% c("g", "m/s2")) {
    errors <- c(errors,
                sprintf("acceleration_units must be 'g' or 'm/s2' (got '%s')",
                        config$accel_units))
  }
  
  # Delimiter
  valid_delims <- c(",", ";", "\t", "tab", "TAB")
  if (!config$csv_delim %in% valid_delims) {
    warnings <- c(warnings,
                  sprintf("Unusual CSV delimiter: '%s' (expected: comma, semicolon, or tab)",
                          config$csv_delim))
  }
  
  # Step frequency range
  if (!is.na(config$sf_freq_min) && !is.na(config$sf_freq_max)) {
    if (config$sf_freq_min >= config$sf_freq_max) {
      errors <- c(errors,
                  sprintf("sf_freq_min (%.1f) must be < sf_freq_max (%.1f)",
                          config$sf_freq_min, config$sf_freq_max))
    }
    if (config$sf_freq_min < 0.3) {
      warnings <- c(warnings,
                    sprintf("sf_freq_min = %.1f Hz is very low (0.3 Hz = 18 steps/min)",
                            config$sf_freq_min))
    }
  }
  
  # ACF max lag sanity check
  if (!is.na(config$acf_max_lag) && !is.na(config$sampling_freq)) {
    acf_lag_seconds <- config$acf_max_lag / config$sampling_freq
    if (acf_lag_seconds < 1.0) {
      warnings <- c(warnings,
                    sprintf("acf_max_lag = %d samples (%.2f s) may be too short to detect stride peaks",
                            config$acf_max_lag, acf_lag_seconds))
    }
    if (acf_lag_seconds > 10.0) {
      warnings <- c(warnings,
                    sprintf("acf_max_lag = %d samples (%.1f s) is unusually long",
                            config$acf_max_lag, acf_lag_seconds))
    }
  }
  
  # AMI parameters
  if (!is.na(config$ami_n_bins) && config$ami_n_bins < 8) {
    warnings <- c(warnings,
                  sprintf("ami_n_bins = %d is very low; 64 is standard for gait",
                          config$ami_n_bins))
  }
  
  # Divergence parameters
  if (!is.na(config$aci_range_end) && !is.na(config$divergence_length) &&
      !is.na(config$samples_per_step)) {
    max_strides_in_curve <- config$divergence_length / (config$samples_per_step * 2)
    if (config$aci_range_end > max_strides_in_curve) {
      errors <- c(errors,
                  sprintf("aci_range_end (%.0f strides) exceeds divergence curve length (%.1f strides). Increase divergence_length or decrease aci_range_end.",
                          config$aci_range_end, max_strides_in_curve))
    }
  }
  
  if (!is.na(config$aci_range_start) && !is.na(config$aci_range_end)) {
    if (config$aci_range_start >= config$aci_range_end) {
      errors <- c(errors,
                  sprintf("aci_range_start (%.0f) must be < aci_range_end (%.0f)",
                          config$aci_range_start, config$aci_range_end))
    }
  }
  
  if (!is.na(config$lds_range_end) && config$lds_range_end <= 0) {
    errors <- c(errors,
                sprintf("lds_range_end must be > 0 (got %.2f)", config$lds_range_end))
  }
  
  # Embedding dimension
  if (!is.na(config$embedding_dim) && (config$embedding_dim < 2 || config$embedding_dim > 15)) {
    warnings <- c(warnings,
                  sprintf("Unusual embedding_dimension: %d (gait standard: 5)",
                          config$embedding_dim))
  }
  
  # Mean period vs samples_per_step consistency
  if (!is.na(config$mean_period) && !is.na(config$samples_per_step)) {
    if (config$mean_period != config$samples_per_step) {
      warnings <- c(warnings,
                    sprintf("mean_period (%d) != samples_per_step (%d). Typically these should match.",
                            config$mean_period, config$samples_per_step))
    }
  }
  
  return(list(
    valid    = length(errors) == 0,
    errors   = errors,
    warnings = warnings
  ))
}


#' Print a formatted summary of configuration parameters
#'
#' @param config Named list of configuration parameters
#'
print_config_summary <- function(config) {
  
  w <- 66  # box width
  rule  <- paste0("  +", paste(rep("-", w), collapse = ""), "+\n")
  hrule <- paste0("  +", paste(rep("=", w), collapse = ""), "+\n")
  line  <- function(...) cat(sprintf("  | %-*s|\n", w - 1, sprintf(...)))
  
  cat("\n")
  cat(hrule)
  line("          RACING Configuration Summary")
  cat(hrule)
  line("Config file:  %s", basename(config$config_file))
  line("Study name:   %s", config$study_name)
  cat(rule)
  line("INPUT DATA")
  line("  Data folder:       %s", config$data_folder)
  line("  Sampling freq:     %d Hz", config$sampling_freq)
  line("  Columns (AP/V/ML): %d / %d / %d",
       config$col_ap, config$col_v, config$col_ml)
  line("  Units: %s  |  Header: %s  |  Delimiter: '%s'",
       config$accel_units,
       ifelse(config$csv_header, "yes", "no"),
       config$csv_delim)
  cat(rule)
  line("PREPROCESSING")
  line("  Tilt correction:   %s", ifelse(config$apply_tilt, "ON", "OFF"))
  line("  Auto orientation:  %s", ifelse(config$auto_orient, "ON", "OFF"))
  cat(rule)
  line("SIGNAL PROCESSING")
  line("  Target:            %d steps x %d samples = %d samples",
       config$target_steps, config$samples_per_step, config$target_samples)
  line("  SF detection:      [%.1f - %.1f] Hz",
       config$sf_freq_min, config$sf_freq_max)
  line("  ACF max lag:       %d samples (%.2f s)",
       config$acf_max_lag, config$acf_max_lag_sec)
  line("  ACF peak dist:     %d samples", config$peak_min_distance)
  cat(rule)
  line("NONLINEAR ANALYSIS")
  line("  AMI:               %d bins, %d lags",
       config$ami_n_bins, config$ami_n_lags)
  line("  Embedding dim:     %d", config$embedding_dim)
  line("  Divergence length: %d samples (%.1f strides)",
       config$divergence_length,
       config$divergence_length / (config$samples_per_step * 2))
  line("  Mean period:       %d samples", config$mean_period)
  line("  LDS range:         0 - %.1f strides (%d steps)",
       config$lds_range_end, config$lds_steps)
  line("  ACI range:         %.0f - %.0f strides",
       config$aci_range_start, config$aci_range_end)
  cat(rule)
  line("OUTPUT")
  line("  Output folder:     %s", config$output_folder)
  line("  Export div curves: %s", ifelse(config$export_div_curves, "YES", "NO"))
  line("  Export resampled:  %s", ifelse(config$export_resampled, "YES", "NO"))
  cat(hrule)
  cat("\n")
}
