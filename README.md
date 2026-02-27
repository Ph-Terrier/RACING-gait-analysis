# RACING: R-based Attractor Complexity INdex for Gait

**R implementation of validated nonlinear gait analysis algorithms for triaxial accelerometer data**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://cran.r-project.org/)
[![CodeOcean](https://img.shields.io/badge/CodeOcean-Reproducible%20Capsule-blue)](https://codeocean.com/capsule/9058247)
[![Zenodo Dataset](https://img.shields.io/badge/Data-Zenodo%2010.5281%2Fzenodo.10148824-blue)](https://doi.org/10.5281/zenodo.10148824)
[![Zenodo Code DOI](https://img.shields.io/badge/Data-Zenodo%2010.5281%2Fzenodo.10148824-blue)](https://doi.org/10.5281/zenodo.18801995)

---

## Overview

RACING is an open-source R implementation of nonlinear gait analysis algorithms originally
developed in MATLAB for the ACIER study (Attractor Complexity Index Empirical Rationalization,
HE-Arc Sante Neuchatel, 2021-2024). It enables researchers to extract validated gait stability
and complexity metrics from triaxial lower-back accelerometer recordings, following FAIR
(Findable, Accessible, Interoperable, Reusable) principles.

The package ports the full ACIER analysis pipeline to R, making these methods accessible to
the international research community without requiring a MATLAB license.

### Key metrics computed

| Metric | Description | Clinical relevance |
|--------|-------------|-------------------|
| **LDS** | Local Divergence Score (short-term, 0-0.5 strides) | Gait stability, fall risk |
| **ACI** | Attractor Complexity Index (long-term, 5-12 strides) | Motor complexity, cognitive engagement |
| **Step/Stride regularity** | ACF-based regularity indices | Gait rhythmicity |
| **Step frequency** | FFT-based cadence detection | Walking speed proxy |
| **RMS ratio** | Mediolateral/norm acceleration ratio | Trunk stability |

### Reference dataset

The ACIER dataset (102 participants, 256 Hz, indoor/outdoor conditions) used to validate
this implementation is publicly available on Zenodo:
[https://doi.org/10.5281/zenodo.10148824](https://doi.org/10.5281/zenodo.10148824)

---

## Requirements

- **R >= 4.0**
- **Operating system:** Windows, macOS, or Linux

### R package dependencies

**Required packages** (the pipeline will not run without these):

```r
install.packages(c(
  "openxlsx",  # Excel output (results file)
  "readxl",    # Excel input (configuration file)
  "gsignal"    # MATLAB-compatible polyphase resampling
))
```

**Strongly recommended** (large performance gains):

```r
install.packages("RANN")  # Fast nearest-neighbor search: ~10x speedup for LDS/ACI
```

**Optional** (used if installed, gracefully skipped otherwise):

```r
install.packages(c(
  "signal",      # Fallback resampling if gsignal is unavailable
  "data.table"   # Faster CSV reading for large datasets
))
```

> **Note on gsignal vs signal:** `gsignal` is strongly recommended because it uses
> polyphase resampling that matches MATLAB's `resample()` function exactly.
> The `signal` package produces slightly different results and is only accepted
> as a fallback. Results may deviate from published ACIER values if `signal` is used.

---

## Installation

Clone this repository:

```bash
git clone https://github.com/Ph-Terrier/RACING-gait-analysis.git
cd RACING-gait-analysis
```

No compilation is required. All scripts are pure R source files.

---

## Quick Start

### Batch processing (recommended)

**Step 1 -- Prepare your configuration**

Copy `config/RACING_config_template.xlsx` to your working folder and rename it
(e.g., `MyStudy_config.xlsx`). Open it in Excel and edit the yellow-highlighted
cells in the `1_Configuration` sheet:

- `data_folder`: path to your CSV files
- `sampling_freq`: your accelerometer sampling rate (Hz)
- `col_ap`, `col_v`, `col_ml`: column indices for each axis (1-based)
- `target_steps`: steps to extract (250 for lab walks)
- `output_folder` and `study_name`: where to save results

**Step 2 -- Run the analysis**

All 13 R scripts must be in the same folder. From the command line (run from
the repository root):

```bash
Rscript R/run_gait_analysis_batch.R MyStudy_config.xlsx
```

Or from RStudio (set your working directory to the `R/` folder first):

```r
setwd("R")
config_path <- "MyStudy_config.xlsx"
source("run_gait_analysis_batch.R")
```

**Step 3 -- Collect results**

The pipeline writes `{study_name}_results.xlsx` to your output folder with four sheets:

- `Metadata`: analysis parameters and run date
- `Results`: one row per participant, 21 metric columns
- `Metrics`: column documentation
- `Notes`: quality control guidance

### Standalone script usage

Individual scripts can be sourced and called independently for specific analyses:

```r
# Step frequency detection only
source("R/step_frequency_detection.R")
result <- step_frequency_fft(acc_vertical, sampling_freq = 256)
cat("Cadence:", result$step_freq_bpm, "steps/min\n")

# Nonlinear metrics only (after preprocessing)
source("R/rosenstein_divergence.R")
source("R/fit_divergence_exponents.R")
div_result <- rosenstein_divergence(signal, m = 5, tau = 15,
                                    max_iter = 1800, mean_period = 75)
exponents <- fit_divergence_exponents(div_result$divergence,
                                      samples_per_step = 75)
cat("LDS:", exponents$LDS, "| ACI:", exponents$ACI, "\n")
```

---

## Input Data Format

RACING expects CSV files with at least three columns of triaxial acceleration data
(anteroposterior, vertical, mediolateral). The axis order and column indices are
specified in the configuration file -- no renaming of files is needed.

**Example CSV structure (default: AP=col 1, V=col 2, ML=col 3):**

```
0.012,-0.987,0.023
0.015,-0.991,0.019
0.011,-0.983,0.027
...
```

Units can be `g` or `m/s2` (set in configuration). The sensor should be placed at
the lower back (L4-L5 level), rigidly attached.

**Minimum signal length:** approximately 350 steps recommended (250 steps extracted
after truncation + buffer for step detection).

---

## Script Architecture

The 13 R scripts are organized in three functional layers:

```
Layer 3: Batch Processing System
  +-- run_gait_analysis_batch.R   Main orchestrator (sources all other scripts)
  +-- read_config.R               Parse Excel configuration
  +-- import_csv_lowback.R        Read CSV files with config
  +-- write_results.R             Save Excel results

Layer 2: Signal Preprocessing Pipeline
  +-- truncate_resample.R         Tilt correction + truncation + resampling
                                  (requires Layer 1 scripts to be sourced first)

Layer 1: Standalone Analysis Functions
  +-- tilt_correction.R           Moe & Nilssen (1998) tilt correction
  +-- step_frequency_detection.R  FFT-based cadence detection
  +-- resample_signal.R           MATLAB-compatible polyphase resampling
  +-- compute_rms.R               RMS norm and ML/norm ratio
  +-- acf_gait_regularity.R       ACF step/stride regularity
  +-- ami.R                       Average Mutual Information (optimal delay tau)
  +-- rosenstein_divergence.R     Rosenstein (1993) divergence curve
  +-- fit_divergence_exponents.R  LDS and ACI extraction from divergence curve
```

Layer 1 scripts have no dependencies on other RACING scripts and can be sourced and
called independently. Layer 2 requires the Layer 1 scripts to be available. The batch
orchestrator (Layer 3) sources all scripts automatically.

### Processing pipeline (per participant)

```
CSV import --> Tilt correction --> Step frequency (FFT)
           --> Signal truncation (250 steps)
           --> Signal resampling (18750 samples = 250 steps x 75 samples/step)
           --> RMS + ACF regularity (on truncated signal)
           --> AMI optimal delay (on resampled signal, 4 axes in parallel)
           --> Rosenstein divergence + LDS/ACI (on resampled signal, 4 axes)
           --> Excel output
```

---

## Key Algorithm Parameters

These parameters are validated against the original MATLAB ACIER implementation
and should not be changed without understanding their scientific basis:

| Parameter (config name) | Default value | Rationale |
|-------------------------|---------------|-----------|
| `target_steps` | 250 | Standardized gait epoch length |
| `samples_per_step` | 75 | 256 Hz / ~3.4 steps/s; determines resampled length |
| `embedding_dim` | 5 | Phase space dimensionality for Rosenstein algorithm |
| `lds_range_end` | 0.5 strides | Short-term divergence window (= 1 step = 75 samples) |
| `aci_range_start` / `aci_range_end` | 5 to 12 strides | Long-term divergence window |
| `ami_n_bins` | 64 | Histogram bins for mutual information estimation |

---

## Reproducibility

To exactly reproduce the ACIER study results:

1. Download the ACIER dataset from Zenodo: [https://doi.org/10.5281/zenodo.10148824](https://doi.org/10.5281/zenodo.10148824)
2. Use the pre-configured `config/RACING_config_template.xlsx` (set `sampling_freq = 256`)
3. Run the batch pipeline on the indoor or outdoor condition folders

A reproducible CodeOcean capsule demonstrating end-to-end processing of the
outdoor older adult condition is available at: [NOT AVAILABLE YET]


---

## Citation


For the underlying ACIER dataset and original MATLAB algorithms, please  cite:

> Piergiovanni, S., & Terrier, P. (2024). Validity of linear and nonlinear measures
> of gait variability to characterize aging gait with a single lower back accelerometer.
> *Sensors*, 24(23), 7427. https://doi.org/10.3390/s24237427

---

## License

This project is licensed under the MIT License -- see [LICENSE](LICENSE) for details.

The original MATLAB algorithms were adapted from scripts available on MATLAB Central
under a BSD license, which is compatible with MIT.

---

## Contact and Support

Philippe Terrier -- HE-Arc Sante, Neuchatel, Switzerland

For bug reports and feature requests, please open an issue on this repository:
[https://github.com/Ph-Terrier/RACING-gait-analysis/issues](https://github.com/Ph-Terrier/RACING-gait-analysis/issues)

---

## Project context

RACING was developed as part of the ORD 2025 grant (Open Research Data, HES-SO)
at HE-Arc Sante Neuchatel. The project objective is to make validated nonlinear gait
analysis methods freely accessible to international researchers following FAIR principles,
enabling reproducibility of the ACIER findings and facilitating new clinical research
on gait stability and fall risk in older adults.
