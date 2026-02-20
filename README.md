R scripts for the article De Benedictis et al. Forest structure and understory functional diversity at multiple scales: The importance of median tree height. Ecological Indicators 183, 114687. https://doi.org/10.1016/j.ecolind.2026.114687

## File description

Run the scripts 1–9 following the order to reproduce the analyses. Script 0 is automatically loaded when necessary.

0. functions: contains functions to calculate diversity indices and additional functions required by other scripts
1. nomenclature: creates lookup table of harmonized species names from transect data
2. data_setup: transect and traits data wrangling
3. resampling: performs computerized resampling of the transect
4. indices: calculates functional diversity indices from data
5. RCindices: calculates indices based on rotated components rather than traits
6. variables: structure and canopy data wrangling, variable selection and clustering
7. models: functional linear models and plots
8. soil_topo: extracts site variables for appendix B
9. complete_case: main analyses done removing missing traits instead of imputing them
10. function_test: checks raoq and qdecomp against known data

Scripts 1-9 can be rendered into reports. Pre-rendered versions are provided in PDF format so that the output can be checked without running them.

## Data requirements

Data can be downloaded at <https://doi.org/10.5281/zenodo.17311232> and put in a "data" folder inside the working directory.

## Zenodo repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15088549.svg)](https://doi.org/10.5281/zenodo.15088549)

