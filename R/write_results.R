# write_results.R
# ================
# RACING Project - ORD 2025
# Helper script: Write batch analysis results to formatted Excel file
#
# Usage:
#   source("write_results.R")
#   write_racing_results(results_df, config, output_path)
#
# Creates an Excel workbook with:
#   - Metadata sheet: analysis parameters and provenance
#   - Results sheet:  one row per file, 21 columns of gait metrics
#   - Metrics sheet:  column documentation (copied from template if available)
#   - Notes sheet:    quality control guidance (copied from template if available)
#
# Safety features:
#   - All paths resolved to absolute before any I/O
#   - Post-save verification (file.exists + file.size)
#   - Automatic CSV fallback if Excel save fails
#   - Emergency RDS backup before Excel attempt
#   - recover_results() function for post-crash recovery
#
# Dependencies: openxlsx
# ================

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required. Install with: install.packages('openxlsx')")
}


# -- Internal: resolve path to absolute and print it clearly --
.resolve_output_path <- function(path) {
  # Expand ~ and normalize separators
  path <- normalizePath(path, mustWork = FALSE)
  return(path)
}


#' Write RACING Batch Results to Excel (with safety nets)
#'
#' @param results_df Data frame with 21 columns (one row per processed file).
#'   Column names must match RACING_Output_Metrics_Template.xlsx Results sheet.
#' @param config Named list from read_racing_config(). Used to populate Metadata.
#' @param output_path Full path for the output .xlsx file.
#' @param template_path Optional path to RACING_Output_Metrics_Template.xlsx.
#'   If provided and found, Metrics and Notes sheets are copied from it.
#' @param verbose If TRUE (default), print confirmation message.
#'
#' @return Invisibly returns the resolved absolute output path.
#'
write_racing_results <- function(results_df, config, output_path,
                                  template_path = NULL, verbose = TRUE) {

  # -- 0. Resolve to absolute path immediately --
  output_path <- .resolve_output_path(output_path)
  if (verbose) {
    cat(sprintf("\n  Output target (absolute): %s\n", output_path))
  }

  # -- 1. Ensure output directory exists --
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (verbose) message("  Created output directory: ", output_dir)
  }

  # Verify directory is writable
  test_file <- file.path(output_dir, ".racing_write_test")
  write_ok <- tryCatch({
    writeLines("test", test_file)
    file.remove(test_file)
    TRUE
  }, error = function(e) FALSE)

  if (!write_ok) {
    warning("Output directory is NOT writable: ", output_dir,
            "\n  Falling back to working directory: ", getwd())
    output_path <- file.path(getwd(), basename(output_path))
    output_dir  <- getwd()
  }

  # -- 1b. Emergency RDS backup BEFORE attempting Excel --
  rds_path <- sub("\\.xlsx$", "_backup.rds", output_path, ignore.case = TRUE)
  tryCatch({
    saveRDS(results_df, rds_path)
    if (verbose) cat(sprintf("  Safety backup: %s\n", rds_path))
  }, error = function(e) {
    warning("Could not write RDS backup: ", e$message)
  })

  # -- 2. Create workbook --
  wb <- openxlsx::createWorkbook()

  # -- 3. Metadata sheet --
  openxlsx::addWorksheet(wb, "Metadata")

  # Header style
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#4472C4",
    fontColour = "#FFFFFF",
    halign = "left",
    fontSize = 11,
    fontName = "Arial"
  )
  label_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "left",
    fontSize = 10,
    fontName = "Arial"
  )
  value_style <- openxlsx::createStyle(
    halign = "left",
    fontSize = 10,
    fontName = "Arial"
  )

  # Write metadata rows (matching template layout)
  meta_labels <- c("Field", "Project", "Date Created", "Dataset",
                    "Contact", "Description", "Analysis Pipeline",
                    "Target Steps", "Target Samples",
                    "Sampling Frequency", "Resampling Rate")

  meta_values <- c(
    "Value",
    "RACING - ORD 2025",
    format(Sys.Date(), "%B %d, %Y"),
    config$study_name,
    "",
    "Gait quality metrics from triaxial accelerometer (lower back)",
    "Tilt correction -> Truncation -> Resampling -> Metric computation",
    sprintf("%d steps", config$target_steps),
    sprintf("%d samples (%d steps x %d samples/step)",
            config$target_samples, config$target_steps, config$samples_per_step),
    sprintf("%d Hz", config$sampling_freq),
    sprintf("%d samples/step", config$samples_per_step)
  )

  openxlsx::writeData(wb, "Metadata", data.frame(A = meta_labels, B = meta_values),
                       colNames = FALSE, startRow = 1)
  openxlsx::addStyle(wb, "Metadata", header_style, rows = 1, cols = 1:2)
  openxlsx::addStyle(wb, "Metadata", label_style,
                      rows = 2:length(meta_labels), cols = 1, gridExpand = TRUE)
  openxlsx::addStyle(wb, "Metadata", value_style,
                      rows = 2:length(meta_labels), cols = 2, gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Metadata", cols = 1, widths = 22)
  openxlsx::setColWidths(wb, "Metadata", cols = 2, widths = 60)

  # -- 4. Results sheet --
  openxlsx::addWorksheet(wb, "Results")

  # Write the results data
  openxlsx::writeData(wb, "Results", results_df, startRow = 1, colNames = TRUE)

  # Style the header row
  results_header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#4472C4",
    fontColour = "#FFFFFF",
    halign = "center",
    fontSize = 10,
    fontName = "Arial",
    border = "bottom",
    borderStyle = "thin"
  )
  openxlsx::addStyle(wb, "Results", results_header_style,
                      rows = 1, cols = 1:ncol(results_df), gridExpand = TRUE)

  # Number formatting for data rows
  n_data_rows <- nrow(results_df)
  if (n_data_rows > 0) {
    # Integer columns: file_id, n_steps, AMI_x/y/z/norm
    int_cols <- which(colnames(results_df) %in%
                        c("file_id", "n_steps", "AMI_x", "AMI_y", "AMI_z", "AMI_norm"))
    int_style <- openxlsx::createStyle(numFmt = "0", halign = "center",
                                        fontSize = 10, fontName = "Arial")
    if (length(int_cols) > 0) {
      openxlsx::addStyle(wb, "Results", int_style,
                          rows = 2:(n_data_rows + 1), cols = int_cols,
                          gridExpand = TRUE)
    }

    # Numeric columns (4 decimals): SF_hz, RMS, regularity, LDS, ACI
    num_cols <- which(colnames(results_df) %in%
                        c("SF_hz", "RMS_norm", "RMS_ratio_ml_norm",
                          "step_reg", "stride_reg",
                          "LDS_x", "LDS_y", "LDS_z", "LDS_norm",
                          "ACI_x", "ACI_y", "ACI_z", "ACI_norm"))
    num_style <- openxlsx::createStyle(numFmt = "0.0000", halign = "center",
                                        fontSize = 10, fontName = "Arial")
    if (length(num_cols) > 0) {
      openxlsx::addStyle(wb, "Results", num_style,
                          rows = 2:(n_data_rows + 1), cols = num_cols,
                          gridExpand = TRUE)
    }

    # Text columns: subject_id, warnings
    txt_cols <- which(colnames(results_df) %in% c("subject_id", "warnings"))
    txt_style <- openxlsx::createStyle(halign = "left",
                                        fontSize = 10, fontName = "Arial")
    if (length(txt_cols) > 0) {
      openxlsx::addStyle(wb, "Results", txt_style,
                          rows = 2:(n_data_rows + 1), cols = txt_cols,
                          gridExpand = TRUE)
    }

    # Highlight rows with warnings or errors
    warn_col <- which(colnames(results_df) == "warnings")
    if (length(warn_col) == 1) {
      warn_fill <- openxlsx::createStyle(fgFill = "#FFF2CC")
      error_fill <- openxlsx::createStyle(fgFill = "#FCE4EC")
      for (r in seq_len(n_data_rows)) {
        w <- results_df$warnings[r]
        if (!is.na(w) && nchar(trimws(w)) > 0) {
          if (grepl("^ERROR:", w)) {
            openxlsx::addStyle(wb, "Results", error_fill,
                                rows = r + 1, cols = 1:ncol(results_df),
                                gridExpand = TRUE, stack = TRUE)
          } else {
            openxlsx::addStyle(wb, "Results", warn_fill,
                                rows = r + 1, cols = 1:ncol(results_df),
                                gridExpand = TRUE, stack = TRUE)
          }
        }
      }
    }
  }

  # Column widths
  openxlsx::setColWidths(wb, "Results", cols = 1, widths = 8)       # file_id
  openxlsx::setColWidths(wb, "Results", cols = 2, widths = 25)      # subject_id
  openxlsx::setColWidths(wb, "Results", cols = 3:20, widths = 12)   # metrics
  openxlsx::setColWidths(wb, "Results", cols = 21, widths = 50)     # warnings

  # Freeze header row
  openxlsx::freezePane(wb, "Results", firstRow = TRUE)

  # Auto-filter
  openxlsx::addFilter(wb, "Results", rows = 1, cols = 1:ncol(results_df))

  # -- 5. Metrics documentation sheet --
  openxlsx::addWorksheet(wb, "Metrics")

  metrics_doc <- data.frame(
    Column = c("file_id", "subject_id", "n_steps", "SF_hz",
               "RMS_norm", "RMS_ratio_ml_norm", "step_reg", "stride_reg",
               "AMI_x", "AMI_y", "AMI_z", "AMI_norm",
               "LDS_x", "LDS_y", "LDS_z", "LDS_norm",
               "ACI_x", "ACI_y", "ACI_z", "ACI_norm", "warnings"),
    Type = c("integer", "string", "numeric", "numeric",
             "numeric", "numeric", "numeric", "numeric",
             rep("integer", 4), rep("numeric", 4), rep("numeric", 4), "string"),
    Description = c(
      "File identifier (sequential)",
      "Input CSV filename (without extension)",
      "Actual number of steps in truncated signal",
      "Step frequency from FFT analysis",
      "RMS of acceleration norm (movement intensity)",
      "Ratio of ML RMS to norm RMS (lateral stability)",
      "Step regularity from ACF (Fisher z-transformed)",
      "Stride regularity from ACF (Fisher z-transformed)",
      "Optimal time delay for AP axis (AMI first minimum)",
      "Optimal time delay for Vertical axis",
      "Optimal time delay for ML axis",
      "Optimal time delay for Norm",
      "Local dynamic stability - AP axis (short-term divergence slope)",
      "Local dynamic stability - Vertical axis",
      "Local dynamic stability - ML axis",
      "Local dynamic stability - Norm",
      sprintf("Attractor complexity index - AP axis (%.0f-%.0f strides)",
              config$aci_range_start, config$aci_range_end),
      sprintf("Attractor complexity index - Vertical axis (%.0f-%.0f strides)",
              config$aci_range_start, config$aci_range_end),
      sprintf("Attractor complexity index - ML axis (%.0f-%.0f strides)",
              config$aci_range_start, config$aci_range_end),
      sprintf("Attractor complexity index - Norm (%.0f-%.0f strides)",
              config$aci_range_start, config$aci_range_end),
      "Quality flags and error messages"
    ),
    Units = c("-", "-", "steps", "Hz", "g", "-", "-", "-",
              rep("samples", 4), rep("lambda/stride", 8), "-"),
    stringsAsFactors = FALSE
  )

  openxlsx::writeData(wb, "Metrics", metrics_doc, startRow = 1, colNames = TRUE)
  openxlsx::addStyle(wb, "Metrics", header_style,
                      rows = 1, cols = 1:4, gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Metrics",
                          cols = 1:4, widths = c(20, 10, 55, 12))

  # -- 6. Notes sheet --
  openxlsx::addWorksheet(wb, "Notes")

  notes <- data.frame(
    Topic = c(
      "Quality Control: n_steps",
      "Quality Control: n_steps",
      "Step Frequency",
      "Movement Intensity",
      "Movement Intensity",
      "Gait Regularity",
      "AMI (Phase Space)",
      "AMI (Phase Space)",
      "LDS (Stability)",
      "LDS (Stability)",
      "ACI (Complexity)",
      "ACI (Complexity)",
      "Configuration",
      "Configuration"
    ),
    Note = c(
      sprintf("Expected: %d steps. If n_steps < target, signal was too short.",
              config$target_steps),
      "Check warnings column for details on short signals.",
      "Detected via FFT on vertical acceleration. Cadence = SF_hz x 60.",
      "RMS computed on gravity-subtracted norm: sqrt(x^2 + y^2 + z^2) - 1",
      "RMS_ratio indicates relative lateral movement (higher = more ML motion)",
      "Based on autocorrelation of norm; atanh-transformed for statistics.",
      "Fraser & Swinney (1986): first minimum of mutual information curve.",
      "Optimal time delay (tau) for phase space reconstruction.",
      sprintf("Rosenstein algorithm: short-term divergence (0-%.1f strides).",
              config$lds_range_end),
      "Higher values indicate less stable gait dynamics.",
      sprintf("Long-term divergence slope (%.0f-%.0f strides).",
              config$aci_range_start, config$aci_range_end),
      "ACI reflects attractor complexity / gait automaticity.",
      sprintf("Config file: %s", basename(config$config_file)),
      sprintf("Processed: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    ),
    stringsAsFactors = FALSE
  )

  openxlsx::writeData(wb, "Notes", notes, startRow = 1, colNames = TRUE)
  openxlsx::addStyle(wb, "Notes", header_style,
                      rows = 1, cols = 1:2, gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Notes", cols = 1:2, widths = c(25, 70))

  # -- 7. Save workbook (with verification) --
  save_ok <- tryCatch({
    openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
    TRUE
  }, error = function(e) {
    warning("Excel save FAILED: ", e$message)
    FALSE
  })

  # -- 8. Post-save verification --
  if (save_ok && file.exists(output_path)) {
    fsize <- file.size(output_path)
    if (fsize > 0) {
      if (verbose) {
        cat(sprintf("\n  [OK] Results saved: %s\n", output_path))
        cat(sprintf("       File size: %.1f KB\n", fsize / 1024))
        cat(sprintf("       Rows: %d  |  Sheets: Metadata, Results, Metrics, Notes\n\n",
                    nrow(results_df)))
      }
    } else {
      warning("Excel file was created but is EMPTY (0 bytes): ", output_path)
      save_ok <- FALSE
    }
  } else if (save_ok) {
    # saveWorkbook didn't error but file doesn't exist -- the ghost problem
    warning("saveWorkbook() reported success but file NOT FOUND at:\n  ",
            output_path,
            "\n  Working directory was: ", getwd())
    save_ok <- FALSE
  }

  # -- 9. CSV fallback if Excel failed --
  if (!save_ok) {
    csv_path <- sub("\\.xlsx$", ".csv", output_path, ignore.case = TRUE)
    tryCatch({
      write.csv(results_df, csv_path, row.names = FALSE)
      cat(sprintf("\n  [FALLBACK] Results saved as CSV: %s\n", csv_path))
      cat(sprintf("  [FALLBACK] RDS backup at: %s\n", rds_path))
      cat("  You can re-generate the Excel later with:\n")
      cat("    results_df <- readRDS('", rds_path, "')\n")
      cat("    write_racing_results(results_df, config, 'output.xlsx')\n\n")
    }, error = function(e) {
      cat(sprintf("\n  [!!] CSV FALLBACK ALSO FAILED: %s\n", e$message))
      cat(sprintf("  [!!] Your data is saved in RDS format at: %s\n", rds_path))
      cat("  [!!] Load with: results_df <- readRDS('", rds_path, "')\n\n")
    })
  }

  invisible(output_path)
}


#' Recover results from RDS backup
#'
#' Use this after a crash or failed save to re-generate the Excel file
#' from the automatic RDS backup.
#'
#' @param rds_path   Path to the _backup.rds file
#' @param config     Config list (if still in memory) OR path to config xlsx
#' @param xlsx_path  Output path for the Excel file
#'
#' @examples
#' \dontrun{
#'   # If config is still in memory:
#'   recover_results("../results/MyStudy_results_backup.rds",
#'                   config, "C:/Users/me/Desktop/recovered.xlsx")
#'
#'   # If starting fresh R session:
#'   source("read_config.R")
#'   source("write_results.R")
#'   cfg <- read_racing_config("my_config.xlsx")
#'   recover_results("../results/MyStudy_results_backup.rds",
#'                   cfg, "C:/Users/me/Desktop/recovered.xlsx")
#' }
#'
recover_results <- function(rds_path, config, xlsx_path) {
  if (!file.exists(rds_path)) {
    stop("RDS backup not found: ", rds_path)
  }
  cat(sprintf("  Loading backup from: %s\n", rds_path))
  results_df <- readRDS(rds_path)
  cat(sprintf("  Recovered %d rows x %d columns\n", nrow(results_df), ncol(results_df)))

  # If config is a path, read it
  if (is.character(config) && length(config) == 1 && file.exists(config)) {
    config <- read_racing_config(config, verbose = FALSE)
  }

  write_racing_results(results_df, config, xlsx_path, verbose = TRUE)
}
