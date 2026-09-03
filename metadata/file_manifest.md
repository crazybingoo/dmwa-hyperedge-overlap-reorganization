# File manifest

This release is organized by manuscript figure, matching the Paper 2 repository pattern.

## Source data

- `Source_Data.xlsx`: aggregate-only source data for Fig. 2–Fig. 7, Supplementary Fig. 1–Supplementary Fig. 4 and Supplementary Table 2, with `README`, `Data_Map` and `Data_Dictionary` worksheets.

## MATLAB

- `matlab/Fig_2/`: macro-state reconstruction.
- `matlab/Fig_3/`: threshold-consistent DMWA and undirected-comparator extraction.
- `matlab/Fig_4/`: directed motif extraction.
- `matlab/Fig_6/`: W-level null controls and S-component ablation.
- `matlab/Supplementary_Fig_1_2/`: reconstruction, synchrony and reciprocity-threshold sensitivity.
- `matlab/Supplementary_Fig_3/`: network-size, replicate and motif-sampling convergence.
- `matlab/core_reference/`: standalone DMWA construction, metric and ablation functions.

## Python

- `python/Fig_3/`: patient-level DMWA comparator summaries.
- `python/Fig_5/`: joint role–reciprocity state-space analysis.
- `python/Fig_6/`: density/order matching, corrected W-level null inference and S-ablation summaries.
- `python/Fig_7/`: simulation-grid fit and ranking stability using the complete-case empirical target.

## R figures

- `r_figures/Fig_2/`–`r_figures/Fig_7/`: main-figure rendering scripts; the Fig. 6 public renderer uses aggregate data only.
- `r_figures/Supplementary_Fig_1_2/`, `r_figures/Supplementary_Fig_3/` and `r_figures/Supplementary_Fig_4/`: supplementary-figure rendering scripts; Supplementary Fig. 4 uses aggregate data only.
- `r_figures/shared/`: shared Paper 1 figure style.

## Excluded by design

Raw/preprocessed SEEG, all row-level human-participant values, clinical patient/seizure tables, identifiers, linkage keys, local paths, generated figures, intermediate CSV/MAT files, Fig. 1 clinical signal data and Supplementary Table 1 are not included.
