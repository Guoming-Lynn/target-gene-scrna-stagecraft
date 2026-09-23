args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
test_root <- normalizePath(file.path(dirname(file_arg), ".."), winslash = "/")
setwd(test_root)
suppressPackageStartupMessages({library(Matrix); library(edgeR); library(limma); library(yaml); library(jsonlite); library(digest)})

root <- tempfile("stagecraft_pathway_")
for (folder in c("00_protocol_manifest", "05_logs", "02_tables", "resources")) {
  dir.create(file.path(root, folder), recursive = TRUE)
}
set.seed(3)
genes <- paste0("G", sprintf("%02d", 1:30))
exposure <- c(0, 1, 2, 3, 4, 5)
counts <- matrix(rnbinom(length(genes) * length(exposure), mu = 30, size = 8), length(genes), length(exposure))
rownames(counts) <- genes
y <- DGEList(counts = counts)
y <- calcNormFactors(y)
design <- model.matrix(~ exposure)
v <- voom(y, design, plot = FALSE)
fit <- eBayes(lmFit(v, design))
saveRDS(list(
  fit = fit, voom = v, design = design, exposure = "exposure",
  consensus = NA_real_, blocked = FALSE, block = paste0("d", seq_along(exposure)),
  y = y, target_gene = "TARGET"
), file.path(root, "05_logs/primary_fit.rds"))

gmt <- setNames(list(genes[1:12], genes[13:24]), c("UP", "DOWN"))
gmt_path <- file.path(root, "resources/toy.json")
write(toJSON(gmt, auto_unbox = TRUE), gmt_path)
cfg <- list(
  random_seed = 3L,
  paths = list(tables_dir = "02_tables", fit_rds = "05_logs/primary_fit.rds"),
  models = list(q_cut = 0.05),
  pathways = list(inter_gene_cor = 0.01, min_size = 10L, max_size = 500L),
  gmt = list(list(path = "resources/toy.json", sha256 = digest(file = gmt_path, algo = "sha256")))
)
cfg_path <- file.path(root, "00_protocol_manifest/analysis_config.yaml")
write_yaml(cfg, cfg_path)
runner <- file.path(test_root, "scripts/part5_run_pathways.R")
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
output <- system2(rscript, c("--vanilla", shQuote(runner), shQuote(cfg_path)), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(output, "status"))) stop(paste(output, collapse = "\n"))
evidence <- read.csv(file.path(root, "02_tables/pathway_evidence.csv"))
audit <- fromJSON(file.path(root, "02_tables/pathway_audit.json"))
stopifnot(
  all(c("q_bh_camera", "q_bh_fgsea", "dual_method_candidate") %in% names(evidence)),
  nrow(evidence) >= 1L,
  isTRUE(all.equal(audit$inter_gene_cor, 0.01)),
  isFALSE(audit$blocked),
  !is.null(audit$packages$fgsea)
)
cat("PATHWAY_CONTRACT_OK\n")
