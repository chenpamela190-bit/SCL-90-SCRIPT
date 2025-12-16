# SCL-90-R analysis scripts

This repository contains analysis scripts to reproduce results reported in:
**[Your manuscript title here]**

## Data
The dataset is publicly available at openICPSR:
https://doi.org/10.3886/E230221V1

## Repository structure
- `scripts/` : analysis scripts
- `outputs/` : generated tables/figures (optional)
- `docs/` : supplementary notes (optional)

## How to reproduce
### Software
- R version: [e.g., 4.4.0]
- Key packages: [e.g., tidyverse, mclust, qgraph, bootnet, NetworkComparisonTest]

### Steps
1. Download data from openICPSR and place it in `data/` (not tracked by Git).
2. Run scripts in this order:
   - `01_qc.R`
   - `02_train_freeze_classifier.R`
   - `03_oos_classification_2022_2024.R`
   - `04_network_estimation_ebicglasso.R`
   - `05_permutation_tests.R`
   - `06_tables_figures.R`
3. Outputs will be written to `outputs/`.

## Notes
- Update file paths in `scripts/config.R` (if provided) before running.
