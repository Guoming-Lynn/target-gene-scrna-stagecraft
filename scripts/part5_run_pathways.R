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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: Rscript part5_run_pathways.R analysis_config.yaml")
cfg_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
cfg <- yaml.load_file(cfg_path)
if (!is.null(cfg$random_seed)) set.seed(cfg$random_seed)
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

tables <- resolve(cfg$paths$tables_dir)
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
outputs <- file.path(tables, c("camera_pathways.csv", "fgsea_pathways.csv", "pathway_evidence.csv"))
if (any(file.exists(outputs))) stop("Refusing to overwrite pathway outputs")

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
  idx[vapply(idx, length, integer(1)) >= 10L & vapply(idx, length, integer(1)) <= 500L]
}

universe <- rownames(fit$coefficients)
camera_rows <- list()
fgsea_rows <- list()
have_fgsea <- requireNamespace("fgsea", quietly = TRUE)

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
  cam <- camera(voom, idx, design, contrast = as.numeric(colnames(design) == exposure), sort = FALSE)
  cam$pathway <- rownames(cam)
  cam$library <- lib
  cam$q_bh <- p.adjust(cam$PValue, method = "BH")
  camera_rows[[length(camera_rows) + 1L]] <- data.frame(
    library = cam$library, pathway = cam$pathway, NGenes = cam$NGenes,
    Direction = cam$Direction, PValue = cam$PValue, q_bh = cam$q_bh,
    status = "SUCCESS", stringsAsFactors = FALSE
  )
  if (have_fgsea) {
    fg <- fgsea::fgseaMultilevel(pathways = gmt[names(idx)], stats = stat, minSize = 10, maxSize = 500)
    fgsea_rows[[length(fgsea_rows) + 1L]] <- data.frame(
      library = lib, pathway = fg$pathway, NES = fg$NES, PValue = fg$pval,
      q_bh = fg$padj, status = "SUCCESS", stringsAsFactors = FALSE
    )
  }
}

camera_tab <- do.call(rbind, camera_rows)
camera_tab$q_bh <- p.adjust(camera_tab$PValue, method = "BH")
write.csv(camera_tab, file.path(tables, "camera_pathways.csv"), row.names = FALSE)

if (length(fgsea_rows)) {
  fgsea_tab <- do.call(rbind, fgsea_rows)
  fgsea_tab$q_bh <- p.adjust(fgsea_tab$PValue, method = "BH")
  write.csv(fgsea_tab, file.path(tables, "fgsea_pathways.csv"), row.names = FALSE)
  merged <- merge(
    camera_tab[, c("library", "pathway", "Direction", "q_bh")],
    fgsea_tab[, c("library", "pathway", "NES", "q_bh")],
    by = c("library", "pathway"), suffixes = c("_camera", "_fgsea")
  )
  merged$same_direction <- (merged$Direction == "Up" & merged$NES > 0) | (merged$Direction == "Down" & merged$NES < 0)
  merged$dual_method_candidate <- merged$q_bh_camera < 0.05 & merged$q_bh_fgsea < 0.05 & merged$same_direction
  merged$robust_primary <- NA
  merged$robustness_status <- "PENDING_PATHWAY_LODO_AND_TECHNICAL_LEADING_EDGE"
  write.csv(merged, file.path(tables, "pathway_evidence.csv"), row.names = FALSE)
  message("Dual-method candidates: ", sum(merged$dual_method_candidate, na.rm = TRUE),
          "; pathway LODO and technical leading-edge audit remain required")
} else {
  write.csv(camera_tab, file.path(tables, "pathway_evidence.csv"), row.names = FALSE)
  message("fgsea package not installed; CAMERA-only table written. Dual-method robust_primary is not available.")
}

message("Empty CAMERA after a large shift is expected. Do not swap the GMT to obtain a paragraph.")

