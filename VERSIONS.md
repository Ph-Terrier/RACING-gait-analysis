# RACING -- Dependency Versions

## Tested environment

The cross-platform validation (MATLAB vs R) reported in the manuscript was
performed with the following package versions:

| Package      | Version | Role                                    | Required |
|--------------|--------:|-----------------------------------------|----------|
| R            |   4.5.1 | Runtime                                 | >= 4.0   |
| gsignal      |  0.3-7  | MATLAB-compatible polyphase resampling  | Yes*     |
| openxlsx     |  4.2.8  | Excel output (results workbook)         | Yes      |
| readxl       |  1.4.5  | Excel input (configuration file)        | Yes      |
| RANN         |  2.6.2  | Fast nearest-neighbor search (KD-tree)  | No**     |
| data.table   | 1.17.6  | Fast CSV reading                        | No**     |

\* Either gsignal or signal is required for resampling. gsignal is strongly
recommended because its polyphase filter matches MATLAB's `resample()`.

\** Optional but recommended for performance. The pipeline falls back to
base R functions when these packages are not installed.

## Version sensitivity

Most RACING dependencies use mature, stable APIs that are robust to routine
CRAN updates:

- **readxl**: only `read_excel()` is used -- stable API since package inception.
- **openxlsx**: standard workbook creation and cell writing -- stable API.
- **RANN**: only `nn2()` is used -- single-purpose package, unchanged API.
- **data.table**: only `fread()` is used -- one of the most stable R functions.

**gsignal is the most version-sensitive dependency.** RACING calls
`gsignal::upfirdn()` and `gsignal::kaiser()` for polyphase resampling that
must produce output numerically equivalent to MATLAB's `resample()`. The
validation in Table 1 of the manuscript was performed with gsignal 0.3-7.
Any future change to the compiled C++ resampling filter in gsignal could
affect numerical results. If in doubt, use the CodeOcean capsule for exact
reproducibility.

## Exact reproducibility

The CodeOcean capsule (ID 9058247) pins all package versions in its
environment and guarantees bitwise-identical results. Use the capsule when
exact reproduction of the published results is required.

For local installations, current CRAN releases are expected to work correctly.

## Full session snapshot

See [`sessionInfo.txt`](sessionInfo.txt) for the complete R session snapshot
(platform, locale, and all loaded package versions).

---

*Last updated: 2026-04-02*
