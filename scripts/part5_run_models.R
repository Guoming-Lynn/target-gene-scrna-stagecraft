#!/usr/bin/env Rscript
# Part 5 gene engine: target-excluded TMM -> voomWithQualityWeights ->
# duplicateCorrelation -> robust eBayes. Support: edgeR QL on FULL.
# Holdout subsets come from analysis_config.yaml. NOT_ESTIMABLE is a result.
#
# Usage:
#   Rscript scripts/part5_run_models.R path/to/00_protocol_manifest/analysis_config.yaml
suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(jsonlite)
  library(yaml)
  library(digest)
})
options(stringsAsFactors = FALSE)
file_arg <- commandArgs()[grepl("^--file=", commandArgs())]
script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg[[1]])) else "scripts"
source(file.path(script_dir, "part5_model_audit.R"))
source(file.path(script_dir, "design_diagnostics.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: Rscript part5_run_models.R analysis_config.yaml")
cfg_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
cfg <- yaml.load_file(cfg_path)
set.seed(cfg$random_seed %||% 42L)
stage_root <- dirname(dirname(cfg_path))

resolve <- function(p) {
  if (is.null(p) || !nzchar(p)) return(NULL)
  if (grepl("^(/|[A-Za-z]:)", p)) return(p)
  file.path(stage_root, p)
}

pb <- resolve(cfg$paths$pseudobulk_dir)
tables <- resolve(cfg$paths$tables_dir)
logs <- resolve(cfg$paths$logs_dir)
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(logs, recursive = TRUE, showWarnings = FALSE)
fit_path <- resolve(cfg$paths$fit_rds %||% "05_logs/primary_fit.rds")
outputs <- c(file.path(tables, c("gene_effects.csv", "gene_evidence.csv", "donor_loo.csv", "fold_audit.csv")),
             file.path(logs, "model_audit.json"), fit_path)
if (any(file.exists(outputs))) stop("Refusing to overwrite existing Part 5 outputs")

target <- as.character(cfg$target_gene)
stopifnot(nzchar(target))
dataset_key <- if (!is.null(cfg$obs$dataset_key)) cfg$obs$dataset_key else "dataset"
block_key <- if (!is.null(cfg$design$block)) cfg$design$block else "dataset_donor_id"
donor_key <- cfg$design$donor_key %||% "dataset_donor_id"
min_rdf <- as.integer(cfg$eligibility$min_rdf_fit)
if (is.na(min_rdf)) min_rdf <- 4L
cpm_cut <- as.numeric(cfg$gene_filter$cpm_threshold)
min_frac <- as.numeric(cfg$gene_filter$min_unit_fraction)
min_abs <- as.integer(cfg$gene_filter$min_units_absolute)
primary_exposure <- if (!is.null(cfg$exposure$primary)) cfg$exposure$primary else "jeffreys_per_10pct"
alt_exposure <- cfg$exposure$alternative
drop_single <- isTRUE(cfg$design$drop_single_level)

meta <- read.csv(file.path(pb, "metadata.csv"), stringsAsFactors = FALSE,
                 check.names = FALSE, fileEncoding = "UTF-8-BOM")
genes <- read.csv(file.path(pb, "genes.csv"), stringsAsFactors = FALSE,
                  check.names = FALSE, fileEncoding = "UTF-8-BOM")
counts <- readMM(file.path(pb, "counts.mtx"))
stopifnot(nrow(counts) == nrow(genes), ncol(counts) == nrow(meta))
if (!"unit_id" %in% names(meta) || anyNA(meta$unit_id) || anyDuplicated(meta$unit_id)) stop("Unique unit_id required")
if (anyNA(genes$gene) || anyDuplicated(genes$gene)) stop("Unique gene names required")
rownames(counts) <- genes$gene
colnames(counts) <- meta$unit_id
if (file.exists(file.path(tables, "metadata_with_source_block.csv"))) {
  annotated <- read.csv(file.path(tables, "metadata_with_source_block.csv"),
                        check.names = FALSE, fileEncoding = "UTF-8-BOM")
  if (!"unit_id" %in% names(annotated) || anyDuplicated(annotated$unit_id)) stop("Unique annotated unit_id required")
  idx <- match(annotated$unit_id, meta$unit_id)
  if (anyNA(idx)) stop("Source metadata contains unknown units")
  original_eligible <- if ("eligible" %in% names(meta)) as.logical(meta$eligible) else rep(TRUE, nrow(meta))
  if (anyNA(original_eligible) || !all(meta$unit_id[original_eligible] %in% annotated$unit_id)) stop("Source metadata omitted eligible units")
  shared <- intersect(names(meta), names(annotated))
  if (!all(vapply(shared, function(k) identical(as.character(meta[[k]][idx]), as.character(annotated[[k]])), logical(1)))) {
    stop("Source metadata changed upstream unit values")
  }
  counts <- counts[, idx, drop = FALSE]
  meta <- annotated
}

if ("eligible" %in% names(meta)) {
  keep_u <- as.logical(meta$eligible)
  if (anyNA(keep_u)) stop("Invalid eligible values")
  meta <- meta[keep_u, , drop = FALSE]
  counts <- counts[, keep_u, drop = FALSE]
}
if (target %in% rownames(counts)) {
  counts <- counts[setdiff(rownames(counts), target), , drop = FALSE]
}

mapping <- yaml.load_file(resolve(cfg$source_block_map))
expected_blocks <- unname(unlist(mapping)[as.character(meta[[dataset_key]])])
if (anyNA(expected_blocks) || any(!nzchar(expected_blocks))) stop("Source map must cover every dataset")
if ("source_block" %in% names(meta) && !identical(as.character(meta$source_block), expected_blocks)) stop("Source map mismatch")
meta$source_block <- expected_blocks
if (length(cfg$arms) != 1L) stop("Run exactly one frozen arm per stage/config")
labels <- unlist(cfg$arms[[1]]$labels)
group_key <- cfg$obs$group_key
keep_arm <- if (length(labels)) meta[[group_key]] %in% labels else
  !grepl("unresolved|stressed|doublet|debris|contaminant|low[ _-]?qc", meta[[group_key]], ignore.case = TRUE)
meta <- meta[keep_arm, , drop = FALSE]
counts <- counts[, keep_arm, drop = FALSE]
if (!nrow(meta)) stop("No eligible units in the frozen arm")
mode <- cfg$estimand$mode
if (is.null(mode) || !mode %in% c("subtype_specific", "joint_common_slope")) stop("Freeze a supported estimand.mode")
for (key in c(donor_key, block_key)) {
  if (!key %in% names(meta) || anyNA(meta[[key]]) || any(!nzchar(as.character(meta[[key]])))) stop("Missing donor/block identity")
}
if (anyDuplicated(meta[, c(donor_key, group_key), drop=FALSE])) stop("Duplicate biological donor-subtype units")
if (mode == "subtype_specific" && (length(unique(meta[[group_key]])) != 1L || anyDuplicated(meta[[donor_key]]))) stop("Subtype-specific mode needs one subtype and one row per donor")
if (anyDuplicated(meta[[donor_key]]) && block_key != donor_key) stop("Repeated donor rows require donor blocking")
if (mode == "joint_common_slope" && !group_key %in% cfg$design$categorical) stop("Joint mode requires subtype fixed effect")
# This runs once on all eligible units in the frozen arm. Every subset inherits it.
meta <- freeze_standardized_covariates(meta, cfg$exposure$freeze_z_on)

placeholder_row <- function(subset, model, status, notes, n_units, n_donors, n_datasets, rdf = NA_integer_) {
  data.frame(
    subset = subset, model = model, gene = NA_character_,
    logFC = NA_real_, CI_L = NA_real_, CI_R = NA_real_, P = NA_real_, q_bh = NA_real_,
    exposure_scale = primary_exposure, n_units = n_units, n_donors = n_donors,
    n_datasets = n_datasets, rdf = rdf, status = status, notes = notes,
    condition_number=NA_real_, exposure_vif=NA_real_, max_exposure_spearman=NA_real_,
    collinearity_review=NA, exposure_sd=NA_real_, exposure_iqr=NA_real_, n_source_blocks=NA_integer_,
    stringsAsFactors = FALSE
  )
}

filter_y <- function(mat) {
  y <- DGEList(counts = mat)
  min_units <- max(min_abs, ceiling(min_frac * ncol(y)))
  keep <- rowSums(cpm(y) >= cpm_cut) >= min_units
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  list(y = y, min_units = min_units)
}

build_design <- function(meta_sub, exposure_name) {
  cat_terms <- cfg$design$categorical
  num_terms <- cfg$design$numeric
  if (is.null(cat_terms)) cat_terms <- character()
  if (is.null(num_terms)) num_terms <- exposure_name
  terms <- c(as.character(cat_terms), as.character(num_terms))
  terms <- unique(c(terms[terms != primary_exposure], exposure_name))
  if (!all(terms %in% names(meta_sub))) stop("Missing frozen design column")
  if (anyNA(meta_sub[, terms, drop = FALSE])) stop("Missing design values")
  for (tm in intersect(terms, as.character(cat_terms))) meta_sub[[tm]] <- factor(meta_sub[[tm]])
  for (tm in setdiff(terms, as.character(cat_terms))) {
    if (!is.numeric(meta_sub[[tm]]) || any(!is.finite(meta_sub[[tm]]))) stop("Non-finite numeric design column")
  }
  if (isTRUE(drop_single)) {
    keep <- vapply(terms, function(tm) {
      if (is.numeric(meta_sub[[tm]])) return(TRUE)
      length(unique(as.character(meta_sub[[tm]]))) >= 2L
    }, logical(1))
    dropped <- terms[!keep]
    terms <- terms[keep]
  } else {
    dropped <- character()
  }
  if (!length(terms)) stop("no design terms remain")
  fml <- as.formula(paste("~", paste(terms, collapse = " + ")))
  design <- model.matrix(fml, data = meta_sub)
  list(design = design, dropped = dropped, terms = terms)
}

fit_blocked <- function(y, design, block) {
  v0 <- voomWithQualityWeights(y, design, plot = FALSE, normalize.method = "none")
  if (length(unique(block)) < 3L || all(table(block) == 1L)) {
    fit <- eBayes(lmFit(v0, design), robust = TRUE)
    return(list(fit = fit, voom = v0, consensus = NA_real_, blocked = FALSE))
  }
  corfit <- duplicateCorrelation(v0, design, block = block)
  v <- voomWithQualityWeights(y, design, block = block, correlation = corfit$consensus,
                              plot = FALSE, normalize.method = "none")
  fit <- eBayes(lmFit(v, design, block = block, correlation = corfit$consensus), robust = TRUE)
  list(fit = fit, voom = v, consensus = corfit$consensus, blocked = TRUE)
}

extract_coef <- function(fit, coef_name, n_units, n_donors, n_datasets, rdf, subset, model, exposure_name, notes) {
  if (!coef_name %in% colnames(fit$coefficients)) {
    return(placeholder_row(subset, model, "NOT_ESTIMABLE",
                           sprintf("coefficient '%s' absent", coef_name),
                           n_units, n_donors, n_datasets, rdf))
  }
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none", confint = TRUE)
  tt$gene <- rownames(tt)
  ci_l <- if ("CI.L" %in% names(tt)) tt$CI.L else tt$logFC
  ci_r <- if ("CI.R" %in% names(tt)) tt$CI.R else tt$logFC
  data.frame(
    subset = subset, model = model, gene = tt$gene,
    logFC = tt$logFC, CI_L = ci_l, CI_R = ci_r,
    P = tt$P.Value, q_bh = tt$adj.P.Val,
    exposure_scale = exposure_name, n_units = n_units, n_donors = n_donors,
    n_datasets = n_datasets, rdf = rdf, status = "SUCCESS", notes = notes,
    condition_number=NA_real_, exposure_vif=NA_real_, max_exposure_spearman=NA_real_,
    collinearity_review=NA, exposure_sd=NA_real_, exposure_iqr=NA_real_, n_source_blocks=NA_integer_,
    stringsAsFactors = FALSE
  )
}

run_limma <- function(meta_sub, counts_sub, exposure_name, subset, model, rdf_floor = min_rdf) {
  n_units <- nrow(meta_sub)
  n_donors <- length(unique(as.character(meta_sub[[donor_key]])))
  n_datasets <- length(unique(as.character(meta_sub[[dataset_key]])))
  sd_floor <- if (exposure_name == "jeffreys_per_10pct") cfg$exposure$sd_jeffreys_floor else if (exposure_name == "log2cpm") cfg$exposure$sd_log2cpm_floor else 1e-8
  exposure_sd <- if (exposure_name %in% names(meta_sub)) sd(meta_sub[[exposure_name]]) else NA_real_
  if (n_donors < (cfg$eligibility$n_exploratory %||% 8L) || !is.finite(exposure_sd) || exposure_sd < (sd_floor %||% 1e-8)) {
    return(list(table = placeholder_row(subset, model, "NOT_ESTIMABLE", "donor/exposure gate failed", n_units, n_donors, n_datasets), fit = NULL))
  }
  built <- tryCatch(build_design(meta_sub, exposure_name), error = function(e) NULL)
  if (is.null(built)) {
    return(list(
      table = placeholder_row(subset, model, "NOT_ESTIMABLE", "design failed", n_units, n_donors, n_datasets),
      fit = NULL
    ))
  }
  design <- built$design
  rank <- qr(design)$rank
  rdf <- nrow(design) - rank
  notes <- if (length(built$dropped)) paste("dropped", paste(built$dropped, collapse = ",")) else ""
  if (rank < ncol(design)) {
    return(list(
      table = placeholder_row(subset, model, "NOT_ESTIMABLE",
                              sprintf("rank %d < %d columns", rank, ncol(design)),
                              n_units, n_donors, n_datasets, rdf),
      fit = NULL
    ))
  }
  if (rdf < rdf_floor) {
    return(list(
      table = placeholder_row(subset, model, "NOT_ESTIMABLE",
                              sprintf("residual df %d < min_rdf %d", rdf, rdf_floor),
                              n_units, n_donors, n_datasets, rdf),
      fit = NULL
    ))
  }
  diagnostic <- diagnose_design(design, exposure_name)
  review <- diagnostic$condition_number > (cfg$diagnostics$condition_review %||% 30) ||
    diagnostic$exposure_vif > (cfg$diagnostics$vif_review %||% 10) ||
    diagnostic$max_exposure_spearman > (cfg$diagnostics$rho_review %||% 0.9)
  if (isTRUE(review)) notes <- paste(notes, "COLLINEARITY_REVIEW_REQUIRED")
  fy <- tryCatch(filter_y(counts_sub), error = function(e) NULL)
  if (is.null(fy) || nrow(fy$y) < 5L) {
    return(list(
      table = placeholder_row(subset, model, "NOT_ESTIMABLE", "too few genes after filter",
                              n_units, n_donors, n_datasets, rdf),
      fit = NULL
    ))
  }
  block <- as.character(meta_sub[[block_key]])
  fitted <- tryCatch(fit_blocked(fy$y, design, block), error = function(e) NULL)
  if (is.null(fitted)) {
    return(list(
      table = placeholder_row(subset, model, "NOT_ESTIMABLE", "limma fit failed",
                              n_units, n_donors, n_datasets, rdf),
      fit = NULL
    ))
  }
  tab <- extract_coef(fitted$fit, exposure_name, n_units, n_donors, n_datasets, rdf,
                      subset, model, exposure_name, notes)
  for (key in names(diagnostic)) tab[[key]] <- diagnostic[[key]]
  tab$collinearity_review <- isTRUE(review)
  tab$exposure_sd <- exposure_sd
  tab$exposure_iqr <- IQR(meta_sub[[exposure_name]])
  tab$n_source_blocks <- length(unique(meta_sub$source_block))
  list(table = tab, fit = fitted, design = design, y = fy$y)
}

run_edger <- function(meta_sub, counts_sub, exposure_name, subset) {
  n_units <- nrow(meta_sub)
  n_donors <- length(unique(as.character(meta_sub[[donor_key]])))
  n_datasets <- length(unique(as.character(meta_sub[[dataset_key]])))
  if (anyDuplicated(meta_sub[[block_key]])) return(placeholder_row(subset, "edgeR_QL", "NOT_ESTIMABLE", "Repeated donor rows require a validated support model", n_units, n_donors, n_datasets))
  built <- tryCatch(build_design(meta_sub, exposure_name), error = function(e) NULL)
  if (is.null(built)) {
    return(placeholder_row(subset, "edgeR_QL", "NOT_ESTIMABLE", "design failed", n_units, n_donors, n_datasets))
  }
  fy <- filter_y(counts_sub)
  y <- estimateDisp(fy$y, built$design)
  fit <- glmQLFit(y, built$design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = exposure_name)
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  data.frame(
    subset = subset, model = "edgeR_QL", gene = rownames(tt),
    logFC = tt$logFC, CI_L = NA_real_, CI_R = NA_real_,
    P = tt$PValue, q_bh = tt$FDR,
    exposure_scale = exposure_name, n_units = n_units, n_donors = n_donors,
    n_datasets = n_datasets, rdf = NA_integer_, status = "SUCCESS", notes = "",
    condition_number=NA_real_, exposure_vif=NA_real_, max_exposure_spearman=NA_real_,
    collinearity_review=NA, exposure_sd=NA_real_, exposure_iqr=NA_real_, n_source_blocks=NA_integer_,
    stringsAsFactors = FALSE
  )
}

subset_units <- function(drop_blocks) {
  if (is.null(drop_blocks) || !length(drop_blocks)) return(rep(TRUE, nrow(meta)))
  !as.character(meta$source_block) %in% as.character(unlist(drop_blocks))
}

subsets <- cfg$subsets
if (is.null(subsets) || !length(subsets)) {
  subsets <- list(FULL = list(drop_blocks = list()))
}
if (!"FULL" %in% names(subsets) || length(subsets$FULL$drop_blocks)) stop("FULL must retain all eligible units")
blocks <- sort(unique(as.character(meta$source_block)))
for (spec in subsets) {
  if (!all(unlist(spec$drop_blocks) %in% blocks)) stop("Unknown source block in holdout config")
}
lodo_names <- setNames(paste0("LODO_", seq_along(blocks)), blocks)
if (any(unname(lodo_names) %in% names(subsets))) stop("LODO_ names are reserved")
if (length(blocks) > 1L) {
  for (b in blocks) subsets[[lodo_names[[b]]]] <- list(drop_blocks = b)
}

all_tabs <- list()
full_fit <- NULL
for (subset_name in names(subsets)) {
  keep <- subset_units(subsets[[subset_name]]$drop_blocks)
  meta_s <- meta[keep, , drop = FALSE]
  counts_s <- counts[, keep, drop = FALSE]
  rdf_floor <- if (subset_name == "FULL") min_rdf else cfg$eligibility$min_rdf_holdout %||% min_rdf
  limma_p <- run_limma(meta_s, counts_s, primary_exposure, subset_name, "limma_primary", rdf_floor)
  all_tabs[[length(all_tabs) + 1L]] <- limma_p$table
  if (identical(subset_name, "FULL") && !is.null(limma_p$fit)) full_fit <- limma_p
  if (!is.null(alt_exposure) && nzchar(alt_exposure) && alt_exposure %in% names(meta_s)) {
    limma_a <- run_limma(meta_s, counts_s, alt_exposure, subset_name, "limma_alternative")
    all_tabs[[length(all_tabs) + 1L]] <- limma_a$table
  }
  if (identical(subset_name, "FULL") && isTRUE(cfg$models$support_edger)) {
    ed <- tryCatch(run_edger(meta_s, counts_s, primary_exposure, subset_name), error = function(e) NULL)
    if (is.null(ed)) ed <- placeholder_row(subset_name, "edgeR_QL", "NOT_ESTIMABLE", "edgeR fit failed", nrow(meta_s), length(unique(meta_s[[donor_key]])), length(unique(meta_s[[dataset_key]])))
    all_tabs[[length(all_tabs) + 1L]] <- ed
  }
}

loo_tabs <- list()
for (donor in unique(as.character(meta[[donor_key]]))) {
  keep <- as.character(meta[[donor_key]]) != donor
  fold <- run_limma(meta[keep, , drop = FALSE], counts[, keep, drop = FALSE],
                    primary_exposure, paste0("LOO_", donor), "limma_primary",
                    cfg$eligibility$min_rdf_holdout %||% min_rdf)
  fold$table$omitted_donor <- donor
  loo_tabs[[length(loo_tabs) + 1L]] <- fold$table
}
write.csv(do.call(rbind, loo_tabs), file.path(tables, "donor_loo.csv"), row.names = FALSE)

for (ds in unique(as.character(meta[[dataset_key]]))) {
  keep <- as.character(meta[[dataset_key]]) == ds
  within <- run_limma(meta[keep, , drop = FALSE], counts[, keep, drop = FALSE],
                      primary_exposure, paste0("WITHIN_", ds), "limma_primary")
  all_tabs[[length(all_tabs) + 1L]] <- within$table
}

gene_table <- do.call(rbind, all_tabs)
write.csv(gene_table, file.path(tables, "gene_effects.csv"), row.names = FALSE)

if (!is.null(full_fit) && !is.null(full_fit$fit)) {
  saveRDS(
    list(
      fit = full_fit$fit$fit, voom = full_fit$fit$voom, design = full_fit$design,
      y = full_fit$y, target_gene = target, exposure = primary_exposure
    ),
    file = fit_path
  )
}

result <- audit_models(gene_table, meta, cfg, subsets, lodo_names, resolve)
audit <- result$audit
write.csv(result$genes, file.path(tables, "gene_evidence.csv"), row.names = FALSE)
write.csv(result$folds, file.path(tables, "fold_audit.csv"), row.names = FALSE)
write(toJSON(audit, auto_unbox = TRUE, pretty = TRUE, null = "null"), file.path(logs, "model_audit.json"))
message("wrote gene_effects.csv and model_audit.json")
message("NOT_ESTIMABLE rows are results. Do not refit a cell Wilcoxon.")

