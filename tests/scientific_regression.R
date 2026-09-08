args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
test_root <- normalizePath(file.path(dirname(file_arg), ".."), winslash = "/")
setwd(test_root)
suppressPackageStartupMessages({library(Matrix); library(yaml); library(jsonlite)})
source(file.path(test_root, "scripts/design_diagnostics.R"))
write_bom_csv <- function(x, path, ...) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write.csv(x, tmp, ...)
  payload <- readBin(tmp, "raw", n = file.info(tmp)$size)
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)
  writeBin(payload, con)
}
set.seed(72)
x <- rnorm(100); z <- rnorm(100)
a <- diagnose_design(model.matrix(~x+z), "x")
b <- diagnose_design(model.matrix(~I(x*1000)+z), "I(x * 1000)")
stopifnot(abs(a$condition_number-b$condition_number)<1e-8, a$exposure_vif < 2)
z <- x + rnorm(100, sd=1e-4)
stopifnot(diagnose_design(model.matrix(~x+z), "x")$exposure_vif > 10)
frozen_z <- freeze_standardized_covariates(
  data.frame(n_cells=c(10,20,100), mean_umi_per_cell=c(1,4,9), log2cpm=c(1,2,8)), "full"
)
stopifnot(isTRUE(all.equal(frozen_z$z_log1p_n_cells,
                            as.numeric(scale(log1p(c(10,20,100)))))))
stopifnot(!isTRUE(all.equal(frozen_z$z_log1p_n_cells[1:2],
                            as.numeric(scale(log1p(c(10,20)))))))
bad_freeze <- try(freeze_standardized_covariates(frozen_z, "subset"), silent=TRUE)
stopifnot(inherits(bad_freeze, "try-error"))

root <- tempfile("stagecraft_regression_")
dir.create(root)
for (p in c("00_protocol_manifest", "03_pseudobulk", "02_tables")) dir.create(file.path(root,p))
cfg <- yaml.load_file("scripts/part5_analysis_config.example.yaml")
stopifnot(is.null(cfg$exposure$collinear_spearman), is.null(cfg$models$q_family))
cfg$arms[[1]]$labels <- "S"
cfg$arms[[1]]$unit_includes_group <- FALSE
cfg$exposure$alternative <- NULL
cfg$design$numeric <- "jeffreys_per_10pct"
cfg$robustness$dominant_block <- "A"
cfg$subsets <- list(FULL=list(drop_blocks=list()))
cfg$paths$pseudobulk_dir <- "03_pseudobulk"
cfg$paths$tables_dir <- "02_tables"
cfg$paths$logs_dir <- "05_logs"
cfg$paths$fit_rds <- "05_logs/primary_fit.rds"
cfg$source_block_map <- "00_protocol_manifest/source_block_map.yaml"
n <- 24L
meta <- data.frame(unit_id=paste0("u",1:n), dataset_donor_id=paste0("d",1:n),
 dataset=rep(c("A","B","C"),each=8), donor_id=paste0("d",1:n), subtype="S",
 jeffreys_per_10pct=rep(seq(0,5,length.out=8),3), eligible=TRUE, n_cells=30)
y <- matrix(rnbinom(120*n,mu=200,size=20),120,n)
y[1:10,] <- matrix(rpois(10*n,lambda=rep(100*exp(meta$jeffreys_per_10pct/3),each=10)),10,n)
writeMM(Matrix(y,sparse=TRUE),file.path(root,"03_pseudobulk/counts.mtx"))
write_bom_csv(data.frame(gene=paste0("G",1:120)),file.path(root,"03_pseudobulk/genes.csv"),row.names=FALSE)
write_bom_csv(meta,file.path(root,"03_pseudobulk/metadata.csv"),row.names=FALSE)
annotated_meta <- transform(meta, source_block = dataset)
write_bom_csv(annotated_meta, file.path(root,"02_tables/metadata_with_source_block.csv"), row.names=FALSE)
write_yaml(list(A="A",B="B",C="C"),file.path(root,"00_protocol_manifest/source_block_map.yaml"))
cfg_path <- file.path(root,"00_protocol_manifest/analysis_config.yaml")
write_yaml(cfg,cfg_path)
runner <- "scripts/part5_run_models.R"
rscript <- file.path(R.home("bin"),"Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"),"Rscript")
output <- system2(rscript,c(shQuote(runner),shQuote(cfg_path)),stdout=TRUE,stderr=TRUE)
cat(paste(output,collapse="\n"),"\n")
stopifnot(is.null(attr(output,"status")))
effects <- read.csv(file.path(root,"02_tables/gene_effects.csv"))
audit <- fromJSON(file.path(root,"05_logs/model_audit.json"))
stopifnot(audit$full_status == "SUCCESS", audit$n_donors == 24L,
 audit$source_evidence == "INTERNAL_SENSITIVITY_ONLY",
 all(c("LODO_1","LODO_2","LODO_3") %in% effects$subset),
 all(is.finite(effects$exposure_vif[effects$subset=="FULL" & effects$model=="limma_primary"])))
cat("SCIENTIFIC_REGRESSION_OK\nArtifacts:",root,"\n")

# Repeat biological donors across subtypes: exercise duplicateCorrelation and
# ensure unblocked edgeR is not accepted as independent support.
joint <- tempfile("stagecraft_joint_regression_")
for (p in c("00_protocol_manifest", "03_pseudobulk", "02_tables"))
  dir.create(file.path(joint,p), recursive=TRUE)
idx <- rep(seq_len(n), each=2L)
joint_meta <- meta[idx,]
joint_meta$unit_id <- paste0("j",seq_len(nrow(joint_meta)))
joint_meta$subtype <- rep(c("TYPE1","TYPE2"),n)
joint_counts <- y[,idx] + matrix(rpois(nrow(y)*length(idx),20),nrow(y),length(idx))
writeMM(Matrix(joint_counts,sparse=TRUE),file.path(joint,"03_pseudobulk/counts.mtx"))
write_bom_csv(data.frame(gene=paste0("G",1:120)),file.path(joint,"03_pseudobulk/genes.csv"),row.names=FALSE)
write_bom_csv(joint_meta,file.path(joint,"03_pseudobulk/metadata.csv"),row.names=FALSE)
write_yaml(list(A="A",B="B",C="C"),file.path(joint,"00_protocol_manifest/source_block_map.yaml"))
cfg$estimand$mode <- "joint_common_slope"
cfg$arms[[1]]$labels <- c("TYPE1","TYPE2")
joint_cfg <- file.path(joint,"00_protocol_manifest/analysis_config.yaml")
write_yaml(cfg,joint_cfg)
output <- system2(rscript,c("--vanilla",shQuote(runner),shQuote(joint_cfg)),stdout=TRUE,stderr=TRUE)
if (!is.null(attr(output,"status"))) stop(paste(output,collapse="\n"))
joint_effects <- read.csv(file.path(joint,"02_tables/gene_effects.csv"))
joint_audit <- fromJSON(file.path(joint,"05_logs/model_audit.json"))
stopifnot(joint_audit$n_donors == n, joint_audit$n_units == 2*n,
          joint_audit$full_status == "SUCCESS",
          all(joint_effects$status[joint_effects$model == "edgeR_QL"] == "NOT_ESTIMABLE"),
          !joint_audit$scientifically_calibrated,
          !joint_audit$formal_gates_pass,
          joint_audit$model_mode_status == "EXPLORATORY_ONLY_CALIBRATION_CONCERN",
          joint_audit$precision_status == "DESCRIPTIVE_CI_PRECISION_ONLY")
cat("JOINT_DONOR_REGRESSION_OK\n")

