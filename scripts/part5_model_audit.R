# Result-derived audit for one frozen Part 5 arm. Shared by runner and tests.
`%||%` <- function(x, fallback) if (is.null(x) || !length(x)) fallback else x

audit_models <- function(tab, meta, cfg, subsets, lodo_names, resolve) {
  full <- tab[tab$subset == "FULL" & tab$model == "limma_primary", , drop = FALSE]
  ok <- nrow(full) > 0L && all(full$status == "SUCCESS")
  evidence <- full[full$status == "SUCCESS" & !is.na(full$gene), , drop = FALSE]
  support <- tab[tab$subset == "FULL" & tab$model == "edgeR_QL" & tab$status == "SUCCESS", , drop = FALSE]
  q_cut <- cfg$models$q_cut %||% 0.05
  floor <- cfg$models$logfc_floor %||% 0.10
  evidence$technical <- grepl("^(MT-|RPS|RPL|HBA[12]$|HBB$|HBD$|HBE1$|HBG[12]$|HBM$|HBQ1$|HBZ$)|^(MALAT1|XIST)$", evidence$gene)
  evidence$q_effect_pass <- is.finite(evidence$q_bh) & evidence$q_bh < q_cut & abs(evidence$logFC) >= floor & !evidence$technical
  support_fc <- support$logFC[match(evidence$gene, support$gene)]
  evidence$dual_same_sign <- is.finite(support_fc) & sign(support_fc) == sign(evidence$logFC)
  blocks <- unique(as.character(meta$source_block))
  fold_names <- if (length(blocks) > 1L) unname(lodo_names[blocks]) else character()
  fold_rows <- lapply(fold_names, function(name) {
    fold <- tab[tab$subset == name & tab$model == "limma_primary", , drop = FALSE]
    data.frame(subset = name, held_out_block = names(lodo_names)[match(name, lodo_names)],
               status = if (nrow(fold) && all(fold$status == "SUCCESS")) "SUCCESS" else "NOT_ESTIMABLE",
               rdf = if (nrow(fold)) fold$rdf[1] else NA_real_)
  })
  folds <- if (length(fold_rows)) do.call(rbind, fold_rows) else data.frame(subset=character(), held_out_block=character(), status=character(), rdf=numeric())
  evidence$lodo_all_same_sign <- rep(length(fold_names) > 0L, nrow(evidence))
  evidence$lodo_n_estimable <- integer(nrow(evidence))
  for (name in fold_names) {
    fold <- tab[tab$subset == name & tab$model == "limma_primary" & tab$status == "SUCCESS", , drop = FALSE]
    fc <- fold$logFC[match(evidence$gene, fold$gene)]
    evidence$lodo_n_estimable <- evidence$lodo_n_estimable + as.integer(is.finite(fc))
    evidence$lodo_all_same_sign <- evidence$lodo_all_same_sign & is.finite(fc) & sign(fc) == sign(evidence$logFC)
  }
  # Fully estimable source holdouts are required for an unrestricted pass.
  evidence$robust_primary <- evidence$q_effect_pass & evidence$dual_same_sign & evidence$lodo_all_same_sign
  candidates <- evidence$q_effect_pass
  primary_candidates <- candidates & evidence$dual_same_sign
  donor_key <- cfg$design$donor_key %||% "dataset_donor_id"
  block_key <- cfg$design$block %||% "dataset_donor_id"
  for (key in c(donor_key, block_key)) {
    if (!key %in% names(meta) || anyNA(meta[[key]]) || any(!nzchar(as.character(meta[[key]])))) stop("Missing donor/block identity")
  }
  n_donors <- length(unique(meta[[donor_key]]))
  n_blocks <- length(unique(meta[[block_key]]))
  dominant <- cfg$robustness$dominant_block
  if (is.null(dominant) && length(subsets$DROP_DOMINANT$drop_blocks) == 1L) dominant <- unlist(subsets$DROP_DOMINANT$drop_blocks)
  if (length(blocks) > 1L && (length(dominant) != 1L || !dominant %in% blocks)) stop("Freeze robustness.dominant_block")
  dominant_fold <- if (length(blocks) > 1L) folds[folds$held_out_block == dominant, , drop = FALSE] else folds[FALSE, ]
  drop_status <- if (nrow(dominant_fold)) dominant_fold$status[1] else "NOT_APPLICABLE"
  source_dependent <- identical(drop_status, "NOT_ESTIMABLE") || length(blocks) < 2L
  if (any(primary_candidates) && nrow(dominant_fold) && drop_status == "SUCCESS") {
    d <- tab[tab$subset == dominant_fold$subset[1] & tab$model == "limma_primary", , drop = FALSE]
    fc <- d$logFC[match(evidence$gene[primary_candidates], d$gene)]
    source_dependent <- any(!is.finite(fc) | sign(fc) != sign(evidence$logFC[primary_candidates]))
  }
  rdf <- if (nrow(full)) full$rdf[1] else NA_real_
  ci_width <- if (all(c("CI_L", "CI_R") %in% names(full))) (full$CI_R - full$CI_L) / 2 else numeric()
  ci_width <- ci_width[is.finite(ci_width) & ci_width >= 0]
  precision <- if (length(ci_width)) list(
    n_genes=length(ci_width), median_ci_half_width=median(ci_width),
    p90_ci_half_width=unname(quantile(ci_width, .9)), effect_floor=floor,
    fraction_ci_half_width_above_effect_floor=mean(ci_width > floor),
    interpretation="Marginal coefficient CI precision; not MDE, achieved power or simultaneous coverage"
  ) else NULL
  mode_exploratory <- identical(cfg$estimand$mode, "joint_common_slope")
  formal <- ok && !mode_exploratory && !source_dependent && nrow(folds) == length(blocks) &&
    nrow(folds) > 0L && all(folds$status == "SUCCESS") && n_donors >= cfg$eligibility$n_formal &&
    length(blocks) >= cfg$eligibility$min_datasets_formal && isTRUE(rdf >= cfg$eligibility$min_rdf_formal)
  reproduction <- cfg$reproduction
  reproduced <- NULL
  reproduction_status <- "NOT_APPLICABLE"
  if (isTRUE(reproduction$required)) {
    parent_path <- resolve(reproduction$path)
    expected <- reproduction$sha256
    if (is.null(expected) || !grepl("^[[:xdigit:]]{64}$", expected) ||
        !file.exists(parent_path) || digest::digest(file=parent_path, algo="sha256") != tolower(expected)) stop("Parent reproduction hash missing/mismatched")
    parent <- read.csv(parent_path, check.names=FALSE)
    if ("subset" %in% names(parent)) parent <- parent[parent$subset == "FULL", , drop=FALSE]
    if ("model" %in% names(parent)) parent <- parent[parent$model == "limma_primary", , drop=FALSE]
    if (!all(c("gene", "logFC") %in% names(parent)) || anyDuplicated(parent$gene)) stop("Invalid parent coefficients")
    anchors <- unlist(reproduction$genes)
    if (!length(anchors)) anchors <- parent$gene
    old <- parent$logFC[match(anchors, parent$gene)]
    new <- evidence$logFC[match(anchors, evidence$gene)]
    reproduced <- length(anchors) > 0L && all(is.finite(old) & is.finite(new) & abs(new-old) <= (reproduction$abs_tol %||% 1e-4))
    reproduction_status <- if (reproduced) "PASS" else "FAILED"
  }
  audit <- list(status="SUCCESS", full_status=if (ok) "SUCCESS" else "NOT_ESTIMABLE",
    calibration_status="NOT_ESTABLISHED_FOR_THIS_ANALYSIS", scientifically_calibrated=FALSE,
    model_mode_status=if (mode_exploratory) "EXPLORATORY_ONLY_CALIBRATION_CONCERN" else "ASSUMPTION_BASED_NOT_CERTIFIED",
    precision_status=if (length(ci_width)) "DESCRIPTIVE_CI_PRECISION_ONLY" else "NOT_ASSESSED",
    precision_summary=precision, protocol_chronology="NOT_VERIFIED",
    external_validation="NOT_ASSESSED", source_evidence="INTERNAL_SENSITIVITY_ONLY",
    estimand_mode=cfg$estimand$mode,
    collinearity_review_required=any(full$collinearity_review %in% TRUE),
    full_reproduced=reproduced, reproduction_status=reproduction_status,
    q_pass=any(candidates), effect_floor_pass=any(candidates),
    dual_method_same_sign=any(candidates) && all(evidence$dual_same_sign[candidates]),
    lodo_all_same_sign=any(primary_candidates) && all(evidence$lodo_all_same_sign[primary_candidates]),
    lodo_complete=nrow(folds) > 0L && all(folds$status == "SUCCESS"),
    source_dependent=source_dependent, drop_dominant_status=drop_status,
    formal_gates_pass=formal, rdf_ok=isTRUE(rdf >= cfg$eligibility$min_rdf_formal),
    n_units=nrow(meta), n_donors=n_donors, n_blocks=n_blocks, donor_key=donor_key, block_key=block_key,
    n_datasets=length(unique(meta[[cfg$obs$dataset_key %||% "dataset"]])), n_source_blocks=length(blocks), rdf=rdf,
    n_robust_primary=sum(evidence$robust_primary),
    evidence_ceiling=cfg$evidence_ceiling, tags=list(), reason="")
  # External diagnostics must provide measured flags; absence stays unknown.
  if (!is.null(cfg$audit_flags)) {
    entry <- cfg$audit_flags
    path <- resolve(entry$path)
    if (is.null(entry$sha256) || !grepl("^[[:xdigit:]]{64}$", entry$sha256) ||
        digest::digest(file=path, algo="sha256") != tolower(entry$sha256)) stop("Diagnostic flags hash mismatch")
    flags <- jsonlite::fromJSON(path, simplifyVector=FALSE)
    allowed <- c("confound_cleared", "sensitivity_pass", "collinear_uninterpretable")
    if (!all(names(flags) %in% allowed) || !all(vapply(flags, function(x) is.logical(x) && length(x)==1L && !is.na(x), logical(1)))) stop("Invalid diagnostic flags")
    audit[names(flags)] <- flags
  }
  list(audit=audit, genes=evidence, folds=folds)
}

