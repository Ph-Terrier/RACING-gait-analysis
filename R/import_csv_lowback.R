# ==============================================================================
# import_csv_lowback.R
# ==============================================================================
# RACING Project - ORD 2025
# Simple CSV import for lower back accelerometer data
#
# Designed to work directly with config from read_racing_config()
# Follows RACING_Batch_Script_Blueprint_v2.md specifications
#
# Usage in batch pipeline:
#   source("read_config.R")
#   source("import_csv_lowback.R")
#   config <- read_racing_config("RACING_config.xlsx")
#   data <- import_csv_lowback("subject_001.csv", config)
#
# Dependencies: data.table (recommended) or base R
# ==============================================================================


#' Import Single CSV Accelerometer File
#'
#' @description
#' Imports a single CSV file containing tri-axial accelerometer data using
#' parameters from RACING configuration Excel file (read via read_config.R).
#' Simple, focused function that integrates directly into the batch pipeline.
#'
#' @param file_path Character string. Full path to the CSV file to import.
#' @param config List. Configuration from read_racing_config(). Must contain:
#'   - csv_delim: delimiter (",", ";", or "\t")
#'   - csv_header: logical, whether file has header row  
#'   - col_ap, col_v, col_ml: column indices (1-based)
#'   - accel_units: "g" or "m/s2"
#'   - sampling_freq: sampling frequency (Hz)
#' @param verbose Logical. Whether to print diagnostic messages (default: FALSE).
#'
#' @return
#' A data.frame with 3 columns (ap, v, ml) in 'g' units, with attributes:
#'   - sampling_rate: Numeric, sampling frequency (Hz)
#'   - units: "g" (always converted to g)
#'   - file: Character, original filename
#'   - n_samples: Integer, number of samples
#'   - import_warnings: Character vector of any warnings
#'
#' @details
#' ## Integration with Batch Pipeline
#'
#' This function is Step 1.A in the RACING batch pipeline:
#' 1. read_racing_config() loads all parameters from Excel
#' 2. import_csv_lowback() reads raw CSV data
#' 3. tilt_correction() preprocesses (if config$apply_tilt == TRUE)
#' 4. truncate_resample() standardizes to target length
#' 5. ... (rest of pipeline)
#'
#' ## Robustness Features
#'
#' - Uses data.table::fread if available (fast, robust encoding handling)
#' - Falls back to read.csv if data.table not installed
#' - Validates column indices before extraction
#' - Converts m/sÂ² to g automatically (dividing by 9.81)
#' - Handles missing values (reports count)
#' - Checks for non-numeric data in acceleration columns
#'
#' ## Error Handling
#'
#' Function stops with clear error message if:
#' - File not found
#' - File is not readable
#' - Specified columns don't exist
#' - Data is not numeric
#'
#' Function warns (but continues) if:
#' - File is suspiciously small (<1 KB)
#' - Missing values detected
#' - Extra columns beyond those specified
#'
#' @examples
#' \dontrun{
#' # Standard usage in batch pipeline
#' config <- read_racing_config("RACING_config.xlsx")
#' data <- import_csv_lowback("001_LB_N_O.csv", config)
#' 
#' # Access data
#' head(data)  # ap, v, ml columns in 'g' units
#' attr(data, "sampling_rate")  # 256
#' attr(data, "n_samples")      # Number of samples
#' }
#'
#' @export
import_csv_lowback <- function(file_path, config, verbose = FALSE) {
  
  # ===========================================================================
  # 1. VALIDATE INPUTS
  # ===========================================================================
  
  if (verbose) message("Importing: ", basename(file_path))
  
  # Check required config parameters
  required_params <- c("csv_delim", "csv_header", "col_ap", "col_v", "col_ml",
                      "accel_units", "sampling_freq")
  missing_params <- setdiff(required_params, names(config))
  if (length(missing_params) > 0) {
    stop("Missing required config parameters: ", 
         paste(missing_params, collapse = ", "),
         "\nEnsure config was created by read_racing_config()")
  }
  
  # Check file exists and is readable
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  if (file.access(file_path, mode = 4) != 0) {
    stop("File is not readable: ", file_path)
  }
  
  # Check file size (warn if suspiciously small)
  file_size <- file.info(file_path)$size
  if (file_size < 1000) {
    warning("File is very small (", file_size, " bytes) - may be empty")
  }
  
  
  # ===========================================================================
  # 2. READ CSV FILE
  # ===========================================================================
  
  warnings_vec <- character(0)
  
  # Handle delimiter (config may have "tab" or "TAB" instead of "\t")
  csv_delim <- config$csv_delim
  if (toupper(csv_delim) == "TAB") {
    csv_delim <- "\t"
  }
  
  # Count comment lines at top of file (lines starting with #)
  # Supports commented CSV headers (e.g. LTMM, any annotated dataset)
  # without requiring an extra config parameter
  first_lines <- readLines(file_path, n = 200L, warn = FALSE)
  n_skip <- sum(cumsum(!grepl("^#", first_lines)) == 0L)
  if (verbose && n_skip > 0L)
    message("  Skipping ", n_skip, " comment line(s) starting with '#'")

  # Try data.table::fread first (fast and handles encoding well)
  if (requireNamespace("data.table", quietly = TRUE)) {
    
    if (verbose) message("  Using data.table::fread for fast import")
    
    raw_data <- tryCatch({
      data.table::fread(
        file = file_path,
        sep = csv_delim,
        header = config$csv_header,
        skip = n_skip,          # skip leading # comment lines
        data.table = FALSE,     # return data.frame
        showProgress = FALSE,
        na.strings = c("", "NA", "NaN", "NULL", "-", "N/A")
      )
    }, error = function(e) {
      stop("Failed to read CSV file: ", file_path, "\n",
           "Error: ", conditionMessage(e), "\n",
           "Check delimiter ('", csv_delim, "') and header setting (",
           config$csv_header, ")")
    })
    
  } else {
    # Fallback to base R read.csv
    if (verbose) message("  Using read.csv (data.table not available)")
    
    raw_data <- tryCatch({
      read.csv(
        file = file_path,
        sep = csv_delim,
        header = config$csv_header,
        comment.char = "#",     # skip # comment lines anywhere in file
        stringsAsFactors = FALSE,
        na.strings = c("", "NA", "NaN", "NULL", "-", "N/A")
      )
    }, error = function(e) {
      stop("Failed to read CSV file: ", file_path, "\n",
           "Error: ", conditionMessage(e))
    })
  }

  if (verbose) {
    message(sprintf("  Read %d rows Ã— %d columns", 
                   nrow(raw_data), ncol(raw_data)))
  }
  
  # Check not empty
  if (nrow(raw_data) == 0) {
    stop("CSV file is empty: ", file_path)
  }
  
  
  # ===========================================================================
  # 3. EXTRACT SPECIFIED COLUMNS
  # ===========================================================================
  
  # Check that specified columns exist
  max_col_needed <- max(config$col_ap, config$col_v, config$col_ml)
  if (max_col_needed > ncol(raw_data)) {
    stop(sprintf(
      "CSV has only %d columns, but config specifies column %d\n",
      ncol(raw_data), max_col_needed),
      "Check column indices in config: AP=%d, V=%d, ML=%d",
      config$col_ap, config$col_v, config$col_ml)
  }
  
  # Extract the 3 specified columns
  data <- data.frame(
    ap = raw_data[, config$col_ap],
    v  = raw_data[, config$col_v],
    ml = raw_data[, config$col_ml],
    stringsAsFactors = FALSE
  )
  
  if (verbose && ncol(raw_data) > 3) {
    message(sprintf("  Extracted columns %d, %d, %d (ignoring %d other columns)",
                   config$col_ap, config$col_v, config$col_ml,
                   ncol(raw_data) - 3))
  }
  
  
  # ===========================================================================
  # 4. VALIDATE AND CONVERT DATA TYPES
  # ===========================================================================
  
  # Ensure all columns are numeric
  for (col_name in c("ap", "v", "ml")) {
    
    if (!is.numeric(data[[col_name]])) {
      
      if (verbose) message("  Converting column '", col_name, "' to numeric")
      
      # Attempt conversion
      original <- data[[col_name]]
      converted <- suppressWarnings(as.numeric(original))
      
      # Check for conversion failures
      n_na_before <- sum(is.na(original))
      n_na_after <- sum(is.na(converted))
      
      if (n_na_after > n_na_before) {
        n_failed <- n_na_after - n_na_before
        warning(sprintf(
          "Column '%s': %d non-numeric value(s) converted to NA",
          col_name, n_failed
        ))
        warnings_vec <- c(warnings_vec,
                         sprintf("%s: %d non-numeric values", col_name, n_failed))
      }
      
      data[[col_name]] <- converted
    }
  }
  
  # Report missing values
  na_counts <- sapply(data, function(x) sum(is.na(x)))
  if (any(na_counts > 0)) {
    pct_na <- 100 * max(na_counts) / nrow(data)
    warning(sprintf(
      "Missing values - AP:%d, V:%d, ML:%d (%.1f%% of data)",
      na_counts[1], na_counts[2], na_counts[3], pct_na
    ))
    warnings_vec <- c(warnings_vec,
                     sprintf("NA values: AP=%d, V=%d, ML=%d", 
                            na_counts[1], na_counts[2], na_counts[3]))
  }
  
  
  # ===========================================================================
  # 5. UNIT CONVERSION (if needed)
  # ===========================================================================
  
  units_out <- "g"  # Always return in 'g' units
  
  if (config$accel_units == "m/s2") {
    
    if (verbose) message("  Converting m/sÂ² to g (Ã· 9.81)")
    
    data$ap <- data$ap / 9.81
    data$v  <- data$v  / 9.81
    data$ml <- data$ml / 9.81
    
  } else if (config$accel_units != "g") {
    # Should not happen if config validation worked
    warning("Unknown acceleration units: '", config$accel_units, 
            "'. Assuming 'g'")
    warnings_vec <- c(warnings_vec, "Unknown units, assumed 'g'")
  }
  
  
  # ===========================================================================
  # 6. ATTACH METADATA AS ATTRIBUTES
  # ===========================================================================
  
  attr(data, "sampling_rate") <- config$sampling_freq
  attr(data, "units") <- units_out
  attr(data, "file") <- basename(file_path)
  attr(data, "n_samples") <- nrow(data)
  attr(data, "import_warnings") <- if (length(warnings_vec) > 0) warnings_vec else NULL
  
  # Store original column indices for traceability
  attr(data, "col_indices") <- c(
    ap = config$col_ap,
    v  = config$col_v,
    ml = config$col_ml
  )
  
  
  # ===========================================================================
  # 7. FINAL REPORT
  # ===========================================================================
  
  if (verbose) {
    duration_sec <- nrow(data) / config$sampling_freq
    message(sprintf("  Import successful: %d samples (%.1f seconds)",
                   nrow(data), duration_sec))
    
    if (length(warnings_vec) > 0) {
      message("  Warnings:")
      for (w in warnings_vec) {
        message("    - ", w)
      }
    }
  }
  
  return(data)
}


# ==============================================================================
# DIAGNOSTIC FUNCTION: Check CSV Structure
# ==============================================================================

#' Check CSV file structure before import
#'
#' @description
#' Diagnostic tool to examine CSV file structure and verify config parameters
#' before attempting import. Useful for troubleshooting import problems.
#'
#' @param file_path Path to CSV file
#' @param config Configuration from read_racing_config()
#' @param n_preview Number of lines to preview (default: 5)
#'
#' @return Invisibly returns a list with diagnostic info
#' @export
check_csv_structure <- function(file_path, config, n_preview = 5) {
  
  cat("\n")
  cat("â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n")
  cat("CSV STRUCTURE CHECK:", basename(file_path), "\n")
  cat("â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n")
  
  # File info
  if (!file.exists(file_path)) {
    cat("âœ— ERROR: File not found!\n")
    return(invisible(NULL))
  }
  
  file_info <- file.info(file_path)
  cat(sprintf("File size: %.1f KB\n", file_info$size / 1024))
  cat(sprintf("Modified:  %s\n\n", file_info$mtime))
  
  # Preview first few lines
  cat("First", n_preview, "lines (raw):\n")
  cat("â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€\n")
  first_lines <- readLines(file_path, n = n_preview, warn = FALSE)
  for (i in seq_along(first_lines)) {
    cat(sprintf("%2d: %s\n", i, first_lines[i]))
  }
  cat("â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€\n\n")
  
  # Try to read with current config
  cat("Attempting import with current config:\n")
  cat(sprintf("  Delimiter: '%s'\n", config$csv_delim))
  cat(sprintf("  Header:    %s\n", config$csv_header))
  cat(sprintf("  Columns:   AP=%d, V=%d, ML=%d\n", 
             config$col_ap, config$col_v, config$col_ml))
  cat(sprintf("  Units:     %s\n\n", config$accel_units))
  
  result <- tryCatch({
    
    data <- import_csv_lowback(file_path, config, verbose = FALSE)
    
    cat("âœ“ Import successful!\n\n")
    cat("Data summary:\n")
    cat(sprintf("  Samples:   %d\n", nrow(data)))
    cat(sprintf("  Duration:  %.1f seconds\n", 
               nrow(data) / config$sampling_freq))
    cat(sprintf("  AP range:  [%.2f, %.2f] g\n", 
               min(data$ap, na.rm=TRUE), max(data$ap, na.rm=TRUE)))
    cat(sprintf("  V range:   [%.2f, %.2f] g\n", 
               min(data$v, na.rm=TRUE), max(data$v, na.rm=TRUE)))
    cat(sprintf("  ML range:  [%.2f, %.2f] g\n", 
               min(data$ml, na.rm=TRUE), max(data$ml, na.rm=TRUE)))
    
    # Check for issues
    cat("\nData quality checks:\n")
    na_count <- sum(is.na(data$ap) | is.na(data$v) | is.na(data$ml))
    if (na_count > 0) {
      cat(sprintf("  âš  Missing values: %d (%.1f%%)\n", 
                 na_count, 100*na_count/nrow(data)))
    } else {
      cat("  âœ“ No missing values\n")
    }
    
    # Check for extreme values (>20g suggests wrong units or sensor saturation)
    extreme <- abs(data$ap) > 20 | abs(data$v) > 20 | abs(data$ml) > 20
    if (any(extreme, na.rm=TRUE)) {
      cat(sprintf("  âš  Extreme values: %d samples >20g (check units?)\n",
                 sum(extreme, na.rm=TRUE)))
    } else {
      cat("  âœ“ All values in plausible range (<20g)\n")
    }
    
    return(invisible(list(
      success = TRUE,
      data = data,
      file_info = file_info
    )))
    
  }, error = function(e) {
    cat("âœ— Import failed!\n")
    cat("Error:", conditionMessage(e), "\n\n")
    
    cat("Troubleshooting suggestions:\n")
    cat("  1. Check delimiter - is it actually '", config$csv_delim, "'?\n")
    cat("  2. Check header setting - does file have column names?\n")
    cat("  3. Check column indices - are they correct?\n")
    cat("  4. Try opening file in text editor to verify structure\n")
    
    return(invisible(list(
      success = FALSE,
      error = conditionMessage(e)
    )))
  })
  
  cat("\nâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n\n")
  
  return(result)
}
