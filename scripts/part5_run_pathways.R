#!/usr/bin/env Rscript
# Part 5 pathway engine: CAMERA on the saved voom design, fgsea on moderated t.
# GMT files are local JSON { "SET_NAME": ["GENE", ...] }; hashes verified here.
#
# Usage:
#   Rscript scripts/part5_run_pathways.R path/to/00_protocol_manifest/analysis_config.yaml
suppressPackageStartupMessages({
  library(limma)
  library(jsonlite)
  library(yaml)
  library(digest)
})
options(stringsAsFactors = FALSE)
`%||%` <- function(x, fallback) if (is.null(x) || !length(x)) fallback else x
technical_pattern <- "^(MT-|RPS|RPL|HBA[12]$|HBB$|HBD$|HBE1$|HBG[12]$|HBM$|HBQ1$|HBZ$)|^(MALAT1|XIST)$"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: Rscript part5_run_pathways.R analysis_config.yaml")
if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("fgsea is required for Part 5 pathways. A CAMERA-only table is not dual-method evidence.")
}
cfg_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
cfg <- yaml.load_file(cfg_path)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
random_seed <- as.integer(cfg$random_seed %||% 42L)
set.seed(random_seed)
stage_root <- dirname(dirname(cfg_path))
resolve <- function(p) {
  if (is.null(p) || !nzchar(p)) return(NULL)
  if (grepl("^(/|[A-Za-z]:)", p)) return(p)
  file.path(stage_root, p)
}

rds_path <- resolve(if (!is.null(cfg$paths$fit_rds)) cfg$paths$fit_rds else file.path("05_logs", "primary_fit.rds"))
if (!file.exists(rds_path)) stop("missing primary_fit.rds; run part5_run_models.R first")
obj <- readRDS(rds_path)
fit <- obj$fit
voom <- obj$voom
design <- obj$design
exposure <- obj$exposure
if (!exposure %in% colnames(design)) stop("exposure coefficient missing from saved design")
pathway_cfg <- cfg$pathways
inter_gene_cor <- pathway_cfg$inter_gene_cor %||% 0.01
fgsea_eps <- pathway_cfg$fgsea_eps %||% 1e-50
fgsea_n_perm_simple <- as.integer(pathway_cfg$fgsea_n_perm_simple %||% 1000L)
min_size <- as.integer(pathway_cfg$min_size %||% 10L)
max_size <- as.integer(pathway_cfg$max_size %||% 500L)
q_cut <- cfg$models$q_cut %||% 0.05
technical_le_max <- 0.50
blocked <- isTRUE(obj$blocked) && is.finite(obj$consensus)

tables <- resolve(cfg$paths$tables_dir)
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
outputs <- file.path(tables, c("camera_pathways.csv", "fgsea_pathways.csv", "pathway_evidence.csv", "pathway_audit.json"))
partials <- sub("(\\.[^.]+)$", ".partial\\1", outputs)
if (any(file.exists(outputs)) || any(file.exists(partials))) stop("Refusing to overwrite pathway outputs")

load_gmt <- function(path) {
  raw <- fromJSON(path, simplifyVector = FALSE)
  if (!is.list(raw) || !length(raw)) stop("GMT JSON must be an object of gene lists: ", path)
  lapply(raw, as.character)
}

gmt_entries <- cfg$gmt
if (is.null(gmt_entries) || !length(gmt_entries)) stop("analysis_config.yaml has no gmt: list")
for (entry in gmt_entries) {
  path <- resolve(entry$path)
  expected <- entry$sha256
  if (is.null(expected) || !grepl("^[[:xdigit:]]{64}$", expected)) stop("Freeze GMT sha256: ", entry$path)
  if (!file.exists(path) || digest(file=path, algo="sha256") != tolower(expected)) stop("GMT hash mismatch: ", entry$path)
}
if (anyDuplicated(vapply(gmt_entries, function(e) basename(e$path), character(1)))) stop("GMT library basenames must be unique")

index_from_gmt <- function(gmt, universe) {
  idx <- lapply(gmt, function(genes) which(universe %in% genes))
  idx[vapply(idx, length, integer(1)) >= min_size & vapply(idx, length, integer(1)) <= max_size]
}

universe <- rownames(fit$coefficients)
camera_rows <- list()
fgsea_rows <- list()

stat <- fit$t[, exposure]
names(stat) <- universe

for (entry in gmt_entries) {
  path <- resolve(entry$path)
  if (is.null(path) || !file.exists(path)) stop("missing GMT: ", entry$path)
  gmt <- load_gmt(path)
  idx <- index_from_gmt(gmt, universe)
  lib <- basename(path)
  if (!length(idx)) {
    camera_rows[[length(camera_rows) + 1L]] <- data.frame(
      library = lib, pathway = NA_character_, NGenes = NA_integer_,
      Direction = NA_character_, PValue = NA_real_, q_bh = NA_real_,
      status = "EMPTY_AFTER_SIZE_FILTER", stringsAsFactors = FALSE
    )
    next
  }
  camera_call <- list(
    voom, idx, design,
    contrast = as.numeric(colnames(design) == exposure),
    sort = FALSE, inter.gene.cor = inter_gene_cor
  )
  if (blocked) {
    camera_call$correlation <- obj$consensus
    camera_call$block <- obj$block
  }
  cam <- do.call(camera, camera_call)
  cam$pathway <- rownames(cam)
  cam$library <- lib
  camera_rows[[length(camera_rows) + 1L]] <- data.frame(
    library = cam$library, pathway = cam$pathway, NGenes = cam$NGenes,
    Direction = cam$Direction, PValue = cam$PValue,
    status = "SUCCESS", inter_gene_cor = inter_gene_cor,
    blocked = blocked, stringsAsFactors = FALSE
  )
  fg <- fgsea::fgseaMultilevel(
    pathways = gmt[names(idx)], stats = stat,
    minSize = min_size, maxSize = max_size,
    eps = fgsea_eps, nPermSimple = fgsea_n_perm_simple
  )
  fgsea_rows[[length(fgsea_rows) + 1L]] <- data.frame(
    library = lib, pathway = fg$pathway, NES = fg$NES, PValue = fg$pval,
    status = "SUCCESS", eps = fgsea_eps, n_perm_simple = fgsea_n_perm_simple,
    leading_edge = vapply(fg$leadingEdge, function(g) paste(g, collapse = ";"), character(1)),
    leading_edge_size = lengths(fg$leadingEdge),
    technical_leading_edge_fraction = vapply(fg$leadingEdge, function(g) if (length(g)) mean(grepl(technical_pattern, g)) else NA_real_, numeric(1)),
    stringsAsFactors = FALSE
  )
}

camera_tab <- do.call(rbind, camera_rows)
camera_tab$q_bh <- p.adjust(camera_tab$PValue, method = "BH")

if (!length(fgsea_rows)) stop("fgsea produced no pathway rows")
fgsea_tab <- do.call(rbind, fgsea_rows)
fgsea_tab$q_bh <- p.adjust(fgsea_tab$PValue, method = "BH")
merged <- merge(
  camera_tab[, c("library", "pathway", "Direction", "q_bh")],
  fgsea_tab[, c("library", "pathway", "NES", "q_bh", "technical_leading_edge_fraction")],
  by = c("library", "pathway"), suffixes = c("_camera", "_fgsea"), sort = FALSE
)
merged <- merged[order(merged$library, merged$pathway, method = "radix"), , drop = FALSE]
rownames(merged) <- NULL
merged$same_direction <- (merged$Direction == "Up" & merged$NES > 0) | (merged$Direction == "Down" & merged$NES < 0)
merged$dual_method_candidate <- merged$q_bh_camera < q_cut & merged$q_bh_fgsea < q_cut & merged$same_direction
merged$technical_leading_edge_pass <- is.finite(merged$technical_leading_edge_fraction) & merged$technical_leading_edge_fraction < technical_le_max
merged$q_cut <- q_cut
merged$robust_primary <- NA
merged$robustness_status <- "PENDING_PATHWAY_LODO"
runtime_packages <- c("limma", "fgsea", "yaml", "jsonlite")
audit_json <- toJSON(list(
  q_cut = q_cut, inter_gene_cor = inter_gene_cor, fgsea_eps = fgsea_eps,
  fgsea_n_perm_simple = fgsea_n_perm_simple, min_size = min_size, max_size = max_size,
  technical_leading_edge_max = technical_le_max,
  blocked = blocked, consensus = obj$consensus %||% NA_real_,
  random_seed = random_seed, rng_kind = RNGkind(),
  n_tied_nonzero_stats = sum(duplicated(stat[stat != 0])),
  R = R.version.string,
  packages = setNames(lapply(runtime_packages, function(p) as.character(packageVersion(p))), runtime_packages)
), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
on.exit(unlink(partials[file.exists(partials)]), add = TRUE)
write.csv(camera_tab, partials[1], row.names = FALSE)
write.csv(fgsea_tab, partials[2], row.names = FALSE)
write.csv(merged, partials[3], row.names = FALSE)
write(audit_json, partials[4])
if (!all(file.rename(partials, outputs))) stop("Could not publish pathway outputs")
message("Dual-method candidates: ", sum(merged$dual_method_candidate, na.rm = TRUE),
        "; pathway LODO and technical leading-edge audit remain required")

message("Empty CAMERA after a large shift is expected. Do not swap the GMT to obtain a paragraph.")

