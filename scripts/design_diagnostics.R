# Diagnostics of the actual fitted matrix, including encoded categorical terms.
diagnose_design <- function(design, exposure_name) {
  keep <- colnames(design) != "(Intercept)"
  x <- design[, keep, drop=FALSE]
  if (any(!is.finite(x)) || !ncol(x)) stop("Invalid design for diagnostics")
  centered <- sweep(x, 2L, colMeans(x))
  norms <- sqrt(colSums(centered^2))
  if (any(norms <= .Machine$double.eps)) {
    return(list(condition_number=Inf, exposure_vif=Inf, max_exposure_spearman=NA_real_))
  }
  z <- sweep(centered, 2L, norms, "/")
  s <- svd(z, nu=0L, nv=0L)$d
  condition <- if (length(s) < ncol(z) || min(s) <= .Machine$double.eps) Inf else max(s)/min(s)
  j <- match(exposure_name, colnames(x))
  if (is.na(j)) stop("Exposure absent from actual design")
  y <- z[, j]
  others <- z[, -j, drop=FALSE]
  if (!ncol(others)) return(list(condition_number=condition, exposure_vif=1, max_exposure_spearman=0))
  rss <- sum(lm.fit(cbind(1, others), y)$residuals^2)
  vif <- if (rss <= .Machine$double.eps) Inf else sum(y^2)/rss
  rho <- max(abs(cor(x[, j], x[, -j, drop=FALSE], method="spearman")))
  list(condition_number=condition, exposure_vif=vif, max_exposure_spearman=rho)
}

# Freeze covariates on the complete eligible frozen arm before any LOO/LODO subset.
freeze_standardized_covariates <- function(meta, freeze_z_on) {
  if (!identical(as.character(freeze_z_on), "full")) {
    stop("exposure.freeze_z_on must be 'full'; subset-specific standardization is forbidden")
  }
  for (item in list(c("z_log1p_n_cells", "n_cells"),
                    c("z_log1p_mean_umi", "mean_umi_per_cell"),
                    c("z_log2cpm", "log2cpm"))) {
    if (item[2] %in% names(meta)) {
      values <- meta[[item[2]]]
      if (item[1] != "z_log2cpm") values <- log1p(values)
      meta[[item[1]]] <- if (isTRUE(sd(values) > 0)) as.numeric(scale(values)) else rep(0, length(values))
    }
  }
  meta
}

