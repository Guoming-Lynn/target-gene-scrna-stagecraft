# R packages for Part 5. Install in R, not pip.
#
# CRAN packages:
# install.packages(c("Matrix", "jsonlite", "yaml", "digest", "statmod"))
#
# Bioconductor packages:
# install.packages("BiocManager")
# BiocManager::install(c("limma", "edgeR", "fgsea"))
#
# `statmod` is required by `calibrate_part5_null.R`; `fgsea` is optional only
# when pathway analysis is not run. `scripts/check_environment.py --stage part5`
# verifies the complete eight-package runtime before a Part 5 run.