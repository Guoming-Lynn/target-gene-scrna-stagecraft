#!/usr/bin/env Rscript
# Execute the actual Part 5 gene runner (including LOO/LODO) on NB global-null data.
# Usage: Rscript --vanilla scripts/calibrate_part5_null.R OUT [REPS=1000]
#        [SEED=271828] [MODE=subtype_specific] [DONORS=18] [SOURCES=3] [GENES=1000]
# This starts at pseudobulk counts. It does not calibrate QC, annotation, pathways
# or Geneformer. No non-null power or MDE is inferred from a null experiment.
suppressPackageStartupMessages({library(Matrix); library(yaml); library(jsonlite); library(digest)})
cli <- commandArgs(trailingOnly=TRUE)
if (!length(cli)) stop("Provide a new output directory")
arg <- function(i, default) if (length(cli) >= i) cli[[i]] else default
out <- path.expand(cli[[1]])
reps <- as.integer(arg(2, 1000)); seed <- as.integer(arg(3, 271828))
mode <- arg(4, "subtype_specific"); donors <- as.integer(arg(5, 18))
sources <- as.integer(arg(6, 3)); ng <- as.integer(arg(7, 1000))
stopifnot(reps >= 1, reps <= 100000, seed >= 1, donors >= 8, sources >= 2,
          donors %% sources == 0, ng >= 100, mode %in% c("subtype_specific", "joint_common_slope"))
filearg <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])
skill <- normalizePath(file.path(dirname(filearg), ".."), winslash="/", mustWork=TRUE)
runner <- file.path(skill, "scripts/part5_run_models.R")
if (dir.exists(out) || file.exists(out)) stop("Output must be new; do not overwrite a calibration")
dir.create(out, recursive=TRUE)
out <- normalizePath(out, winslash="/", mustWork=TRUE)
hash_files <- c("scripts/calibrate_part5_null.R", "scripts/part5_run_models.R",
                "scripts/part5_model_audit.R", "scripts/design_diagnostics.R",
                "scripts/part5_analysis_config.example.yaml")
manifest <- list(status="FROZEN_BEFORE_EXECUTION", scenario="GLOBAL_NULL_NB",
  replicates=reps, base_seed=seed, estimand=mode, donors=donors, sources=sources,
  genes=ng, truth_log2_slope=0, q_cut=.05, raw_p_cut=.05, interval_level=.95,
  dispersion=.15, donor_log_sd=.35, source_log_sd=.3, exposure_source_shift=.4,
  generator="Balanced sources; source-shifted exposure independent of gene-specific donor effects; NB counts; repeated subtypes share donor intercepts",
  scope="Pseudobulk gene engine + source LODO + donor LOO; excludes upstream QC/selection, pathways and Part 6",
  code_sha256=setNames(lapply(hash_files, function(p) digest(file=file.path(skill,p), algo="sha256")), hash_files),
  runtime=list(R=R.version.string, packages=setNames(lapply(c("limma","edgeR","Matrix","statmod"),
    function(p) as.character(packageVersion(p))), c("limma","edgeR","Matrix","statmod"))))
write_json(manifest, file.path(out,"manifest.json"), auto_unbox=TRUE, pretty=TRUE)
records <- list()
for (iteration in seq_len(reps)) {
  # Each draw is independent of RNG consumed by the pipeline.
  set.seed(seed + iteration - 1L)
  stage <- file.path(out, sprintf("rep_%05d", iteration))
  for (folder in c("00_protocol_manifest", "03_pseudobulk", "02_tables", "05_logs"))
    dir.create(file.path(stage,folder), recursive=TRUE)
  source_id <- rep(seq_len(sources), each=donors/sources)
  x <- as.numeric(scale(rnorm(donors) + .4 * source_id))
  subtype_n <- if (mode == "joint_common_slope") 2L else 1L
  donor_index <- rep(seq_len(donors), each=subtype_n)
  n <- length(donor_index)
  meta <- data.frame(unit_id=paste0("u",seq_len(n)),
    dataset_donor_id=paste0("d",donor_index), donor_id=paste0("d",donor_index),
    dataset=paste0("SRC", source_id[donor_index]),
    subtype=rep(paste0("TYPE",seq_len(subtype_n)), donors),
    exposure=x[donor_index], eligible=TRUE, n_cells=40)
  baseline <- runif(ng, log(50), log(500))
  donor_effect <- matrix(rnorm(ng*donors, sd=.35), ng, donors)
  source_effect <- matrix(rnorm(ng*sources, sd=.3), ng, sources)
  # No exposure term: the conditional gene slope is zero in both modes.
  log_mu <- baseline + donor_effect[,donor_index] + source_effect[,source_id[donor_index]]
  log_mu <- sweep(log_mu, 2, rnorm(n, sd=.2), "+")
  counts <- matrix(rnbinom(ng*n, mu=as.vector(exp(log_mu)), size=1/.15), ng, n)
  writeMM(Matrix(counts, sparse=TRUE), file.path(stage,"03_pseudobulk/counts.mtx"))
  write.csv(data.frame(gene=paste0("FEATURE",seq_len(ng))), file.path(stage,"03_pseudobulk/genes.csv"), row.names=FALSE)
  write.csv(meta, file.path(stage,"03_pseudobulk/metadata.csv"), row.names=FALSE)
  cfg <- yaml.load_file(file.path(skill,"scripts/part5_analysis_config.example.yaml"))
  cfg$target_gene <- "TARGET_GENE"; cfg$random_seed <- seed + iteration - 1L
  cfg$estimand$mode <- mode
  cfg$arms[[1]]$labels <- unique(meta$subtype)
  cfg$exposure$primary <- "exposure"; cfg$exposure$alternative <- NULL
  cfg$design$numeric <- "exposure"
  cfg$design$categorical <- c("dataset", if (subtype_n == 2L) "subtype")
  cfg$robustness$dominant_block <- "SRC1"
  cfg$subsets <- list(FULL=list(drop_blocks=list()))
  cfg$evidence_ceiling <- "exploratory"
  cfg_path <- file.path(stage,"00_protocol_manifest/analysis_config.yaml")
  write_yaml(cfg, cfg_path)
  ids <- paste0("SRC",seq_len(sources))
  write_yaml(setNames(as.list(ids),ids), file.path(stage,cfg$source_block_map))
  engine <- new.env(parent=globalenv())
  # Same source file and arguments as the CLI; loaded once per fresh environment
  # to avoid repeatedly starting R. Model, filtering and audit code are unmodified.
  engine$commandArgs <- function(trailingOnly=FALSE) {
    if (trailingOnly) cfg_path else c(paste0("--file=",runner),cfg_path)
  }
  error_message <- ""
  logcon <- file(file.path(stage,"05_logs/runner.log"), open="wt")
  sink(logcon, type="message")
  start <- proc.time()[[3]]
  completed <- tryCatch({sys.source(runner,envir=engine); TRUE},
                       error=function(e) {error_message <<- conditionMessage(e); FALSE})
  sink(type="message"); close(logcon)
  row <- data.frame(replicate=iteration, seed=seed+iteration-1L, engine_completed=completed,
    primary_estimable=FALSE, n_tested=0L, raw_fpr=NA_real_, bh_fdp=0,
    ci_coverage=NA_real_, robust_any=FALSE, lodo_failure=TRUE,
    formal_eligibility=FALSE, seconds=proc.time()[[3]]-start, error=error_message)
  if (completed) {
    tab <- read.csv(file.path(stage,"02_tables/gene_effects.csv"))
    full <- tab[tab$subset == "FULL" & tab$model == "limma_primary" & tab$status == "SUCCESS",]
    row$primary_estimable <- nrow(full) > 0
    row$n_tested <- nrow(full)
    if (nrow(full)) {
      row$raw_fpr <- mean(full$P < .05)
      # All-null FDR equals Pr(any BH rejection), not the fraction of rejected genes.
      row$bh_fdp <- as.numeric(any(full$q_bh < .05))
      row$ci_coverage <- mean(full$CI_L <= 0 & full$CI_R >= 0)
    }
    audit <- fromJSON(file.path(stage,"05_logs/model_audit.json"))
    row$lodo_failure <- !isTRUE(audit$lodo_complete)
    row$formal_eligibility <- isTRUE(audit$formal_gates_pass)
    row$robust_any <- audit$n_robust_primary > 0
  }
  records[[iteration]] <- row
  write.csv(do.call(rbind,records), file.path(out,"replicates.csv"), row.names=FALSE)
  if (iteration %% 10L == 0L || iteration == 1L) cat("Completed",iteration,"/",reps,"\n")
}
results <- do.call(rbind,records)
summary_metric <- function(x, bernoulli=FALSE) {
  x <- x[is.finite(x)]; n <- length(x)
  if (!n) return(list(estimate=NULL, n_replicates=0))
  ci <- if (bernoulli) unname(binom.test(sum(x), n)$conf.int) else NULL
  list(estimate=mean(x), mcse=if(n>1) sd(x)/sqrt(n) else NULL,
       binomial_95_ci=ci, n_replicates=n)
}
report <- list(status=if(all(results$engine_completed)) "COMPLETED_LIMITED_NULL_EXPERIMENT" else "COMPLETED_WITH_FAILURES",
  scientifically_calibrated=FALSE, scenario="GLOBAL_NULL_NB", mode=mode,
  requested_replicates=reps, completed_replicates=sum(results$engine_completed),
  raw_type_I_error=summary_metric(results$raw_fpr),
  BH_global_null_FDR=summary_metric(results$bh_fdp,TRUE),
  marginal_CI_coverage=summary_metric(results$ci_coverage),
  NOT_ESTIMABLE_rate=summary_metric(as.numeric(!results$primary_estimable),TRUE),
  LODO_failure_rate=summary_metric(as.numeric(results$lodo_failure),TRUE),
  robust_any_null_rate=summary_metric(as.numeric(results$robust_any),TRUE),
  power=NULL, MDE=NULL,
  limitations=c("Null experiment only; no non-null power or MDE",
    "Repeated gene statistics are dependent; MCSE uses independent replicate means",
    "FDR includes failed fits as no discovery; failure rate is reported separately",
    "Raw FPR/coverage conditional on estimable filtered genes; their effective repetitions are reported",
    "Does not validate other donor/source sizes, depth collinearity, attrition, QC/selection, pathways or Part 6"))
write_json(report, file.path(out,"calibration_report.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=8)
cat("Report:",file.path(out,"calibration_report.json"),"\n")
if (!all(results$engine_completed)) quit(status=2)
