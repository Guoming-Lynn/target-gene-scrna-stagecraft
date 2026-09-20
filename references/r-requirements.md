# R packages for Part 5. Install in R, not pip.
#
# CRAN packages:
# install.packages(c("Matrix", "jsonlite", "yaml", "digest", "statmod"))
#
# Bioconductor packages:
# install.packages("BiocManager")
# BiocManager::install(c("limma", "edgeR", "fgsea"))
#
# `fgsea` is required by `scripts/check_environment.py --stage part5` and by
# `scripts/part5_run_pathways.R`. Do not treat it as optional for a Part 5 run.