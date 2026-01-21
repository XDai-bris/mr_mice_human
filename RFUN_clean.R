# RFUN_clean.R
# ----------------------------------------------------------------------------
# Design goals
# — Tidy, modular utilities for your DO/LOCO/QTL & MR workflows
# - Mouse→Human MR API kept (OpenGWAS IDs input + optional labels)
# ----------------------------------------------------------------------------

# ---------------------------- Dependencies ----------------------------------
require_pkgs <- function(pkgs = c(
  "dplyr","data.table","tibble","stringr","purrr","ggplot2","patchwork",
  "MASS","lme4qtl","lme4","lattice","Matrix","tidyr","pbapply","ggrepel",
  "viridis","qtl2","TwoSampleMR","ieugwasr","rlang","parallel"
)) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop(
    "Missing packages: ", paste(miss, collapse = ", "),
    ". Install with install.packages() or remotes::install_github() as needed."
  )
  invisible(TRUE)
}

# Call once at load time
require_pkgs()

# ---------------------------- Small utilities -------------------------------
# Inverse normal transform (handles NAs)
inverse_normal_transform <- function(x, offset = 0.5, ties.method = "average") {
  stopifnot(is.numeric(x))
  not_na <- is.finite(x)
  r <- rank(x[not_na], ties.method = ties.method)
  n <- sum(not_na)
  u <- (r - offset) / (n - 2 * offset + 1)
  z <- qnorm(u)
  out <- rep(NA_real_, length(x)); out[not_na] <- z
  out
}

# Pick contiguous region around a peak above LOD threshold
.get_one_region <- function(a, index, lod_thresh = 4) {
  r1 <- which(a[index:1] < lod_thresh)[1] - 2
  r2 <- which(a[index:length(a)] < lod_thresh)[1] - 2
  a[(index - r1):(index + r2)] <- 0
  list(a = a, index = index)
}

# Identify all independent tophits indices >= lod_thresh
get_tophits <- function(a, lod_thresh = 4) {
  a[a < lod_thresh] <- 0
  peaks <- integer()
  while (any(a >= lod_thresh)) {
    index <- which.max(a)
    res <- .get_one_region(a, index, lod_thresh)
    peaks <- c(peaks, res$index)
    a <- res$a
  }
  peaks
}

# For one trait in a named list of LOD scans, return marker IDs exceeding threshold
get_markers_for_trait <- function(trait, scan_list, lod_thresh = 4, lod_col = 1) {
  out <- scan_list[[trait]]
  if (is.null(out) || nrow(out) == 0) return(character())
  # choose column
  lod <- if (is.character(lod_col)) out[, lod_col] else out[, lod_col]
  idx <- tryCatch(get_tophits(lod, lod_thresh), error = function(e) integer(0))
  if (!length(idx)) character() else rownames(out)[idx]
}

# Fast univariate association (OLS closed-form)
fast_assoc <- function(y, x) {
  valid <- is.finite(y) & is.finite(x)
  y <- y[valid]; x <- x[valid]
  n <- length(y)
  vx <- var(x); vy <- var(y)
  bhat <- cov(y, x) / vx
  ahat <- mean(y) - bhat * mean(x)
  rsq <- (bhat * vx)^2 / (vx * vy)
  fval <- rsq * (n - 2) / (1 - rsq)
  tval <- sqrt(fval)
  se <- abs(bhat / tval)
  pval <- pf(fval, 1, n - 2, lower.tail = FALSE)
  tibble::tibble(ahat = ahat, bhat = bhat, se = se, fval = fval, pval = pval, n = n)
}

# Vectorised effect extraction across markers (parallelisable)
get_effects_for_trait <- function(g, phen, markers, trait, mc.cores = 1L) {
  idx <- match(markers, rownames(g))
  if (anyNA(idx)) {
    warning("Markers missing in genotype matrix: ", paste(markers[is.na(idx)], collapse = ", "))
    markers <- markers[!is.na(idx)]; idx <- idx[!is.na(idx)]
  }
  if (!length(idx)) return(dplyr::tibble())
  g_inst <- g[idx, , drop = FALSE]
  parallel::mclapply(seq_len(nrow(g_inst)), function(i) {
    fast_assoc(phen[[trait]], g_inst[i, ]) %>%
      dplyr::mutate(marker = rownames(g_inst)[i], trait = trait)
  }, mc.cores = mc.cores) |>
    dplyr::bind_rows()
}

# LOD at a given chr/pos
get_lod_at_pos <- function(scan_obj, map_list, chr, pos, tol = 1e-4) {
  if (!chr %in% names(map_list)) return(NA_real_)
  map_chr <- map_list[[chr]]
  mrk <- names(map_chr)[abs(map_chr - pos) < tol]
  if (!length(mrk)) return(NA_real_)
  as.numeric(scan_obj[mrk, 1])
}

# ------------------------ LOCO LMM (slow/accurate) ---------------------------
# Per-marker LOCO LMM with robust sanity checks
run_loco_lmm_for_markers <- function(
    trait, markers, G, phen, covar_mat, K_loco, snp_chr,
    maf_min = 0.01, maxfun = 2e5, nugget = 1e-6, verbose = TRUE
) {
  N <- length(markers)
  out <- vector("list", N)
  
  make_na_row <- function(s, chr, note, maf = NA_real_) {
    data.frame(snp = s, chr = chr, maf = maf,
               beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
               beta_orig = NA_real_, se_orig = NA_real_, note = note,
               stringsAsFactors = FALSE)
  }
  
  for (i in seq_len(N)) {
    s <- markers[i]; chr <- as.character(snp_chr[s])
    if (verbose) cat(sprintf("[%d/%d] %s (chr %s)\n", i, N, s, ifelse(is.na(chr) || !nzchar(chr), "NA", chr)))
    
    if (!(s %in% colnames(G))) { out[[i]] <- make_na_row(s, chr, "snp_not_in_G"); next }
    if (is.na(chr) || !nzchar(chr)) { out[[i]] <- make_na_row(s, NA, "no_chr_mapping"); next }
    Kc <- K_loco[[chr]]
    if (is.null(Kc)) { out[[i]] <- make_na_row(s, chr, "no_K_for_chr"); next }
    
    ids <- rownames(Kc)
    y  <- phen[[trait]][match(ids, phen$Mouse.ID)]
    X  <- as.data.frame(covar_mat[ids, , drop = FALSE])
    X  <- X[, setdiff(colnames(X), trait), drop = FALSE]
    
    keep <- stats::complete.cases(y, X); y <- y[keep]; X <- X[keep, , drop = FALSE]; Kc <- Kc[keep, keep]
    
    if (ncol(X)) {
      num <- vapply(X, is.numeric, logical(1))
      if (any(num)) X[, num] <- lapply(X[, num, drop = FALSE], function(x) as.numeric(scale(x, TRUE, TRUE)))
    }
    
    g0 <- as.numeric(G[rownames(Kc), s])
    if (all(is.na(g0))) { out[[i]] <- make_na_row(s, chr, "g_all_NA"); next }
    keep_snp <- stats::complete.cases(g0); g0 <- g0[keep_snp]; y <- y[keep_snp]; X <- X[keep_snp,,drop=FALSE]; Kc <- Kc[keep_snp, keep_snp]
    
    m <- mean(g0, na.rm = TRUE); maf <- min(m/2, 1 - m/2)
    if (is.na(maf) || maf < maf_min || stats::sd(g0) == 0) { out[[i]] <- make_na_row(s, chr, "low_var_or_maf", maf); next }
    
    y_scaled <- scale(y, TRUE, TRUE)[,1]; sd_y <- stats::sd(y, na.rm = TRUE)
    if (ncol(X)) {
      num <- vapply(X, is.numeric, logical(1)); if (any(num)) X[, num] <- lapply(X[, num, drop=FALSE], function(x) as.numeric(scale(x, TRUE, TRUE)))
    }
    g_c <- scale(g0, TRUE, FALSE)[,1]
    
    df <- cbind.data.frame(ID = factor(rownames(Kc), levels = rownames(Kc)), y = y_scaled, X, g = g_c)
    drop1 <- vapply(df, function(z) is.factor(z) && nlevels(z) < 2, logical(1)); if (any(drop1)) df <- df[, !drop1, drop = FALSE]
    
    eigmin <- tryCatch(min(eigen((Kc + t(Kc))/2, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA_real_)
    if (!is.na(eigmin) && eigmin < 1e-8) diag(Kc) <- diag(Kc) + nugget
    
    ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = maxfun))
    fit <- tryCatch(lme4qtl::relmatLmer(y ~ g + . - y + (1|ID), data = df, relmat = list(ID = Kc), control = ctrl), error = function(e) e)
    if (inherits(fit, "error")) { out[[i]] <- make_na_row(s, chr, paste0("fit_error: ", fit$message), maf); next }
    
    be <- tryCatch(lme4::fixef(fit)["g"], error = function(e) NA_real_)
    se <- tryCatch(sqrt(diag(vcov(fit)))["g"], error = function(e) NA_real_)
    z  <- if (is.finite(be) && is.finite(se) && se > 0) as.numeric(be/se) else NA_real_
    p  <- if (is.finite(z)) 2*stats::pnorm(-abs(z)) else NA_real_
    
    beta_orig <- be * sd_y; se_orig <- se * sd_y
    
    out[[i]] <- data.frame(snp = s, chr = chr, maf = maf,
                           beta = be, se = se, z = z, p = p,
                           beta_orig = beta_orig, se_orig = se_orig,
                           note = NA_character_, stringsAsFactors = FALSE)
  }
  
  ans <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  if (!is.null(ans$p)) ans$q <- stats::p.adjust(ans$p, method = "BH")
  ans
}

# Fast batched LOCO-LMM by chromosome (re-implemented for speed and stability)
run_fast_loco_lmm <- function(
    trait, markers, G, phen, covar_mat, K_loco, snp_chr,
    maf_min = 0.01, maxfun = 2e5, nugget = 1e-6, verbose = TRUE
) {
  N <- length(markers)
  res <- vector("list", N)
  chr_vec <- as.character(snp_chr[markers])
  idx_by_chr <- split(seq_along(markers), chr_vec)
  
  .i <- 0L
  for (cc in names(idx_by_chr)) {
    snp_idx <- idx_by_chr[[cc]]
    if (is.na(cc) || !nzchar(cc)) {
      for (j in snp_idx) { s <- markers[j]; .i <- .i + 1L; if (verbose) cat(sprintf("[%d/%d] %s (chr NA)\n", .i, N, s))
      res[[j]] <- data.frame(snp=s, chr=NA, maf=NA, beta=NA, se=NA, z=NA, p=NA, beta_orig=NA, se_orig=NA, q=NA, note="no_chr_mapping") }
      next
    }
    if (verbose) cat("\n--- Chromosome ", cc, "---\n", sep = "")
    Kc <- K_loco[[cc]]
    if (is.null(Kc)) {
      for (j in snp_idx) { s <- markers[j]; .i <- .i + 1L; if (verbose) cat(sprintf("[%d/%d] %s (chr %s)\n", .i, N, s, cc))
      res[[j]] <- data.frame(snp=s, chr=cc, maf=NA, beta=NA, se=NA, z=NA, p=NA, beta_orig=NA, se_orig=NA, q=NA, note="no_K_for_chr") }
      next
    }
    
    idsK <- rownames(Kc)
    ids_common <- Reduce(intersect, list(idsK, phen$Mouse.ID, rownames(covar_mat)))
    Kc  <- Kc[ids_common, ids_common, drop = FALSE]
    y0  <- phen[[trait]][match(ids_common, phen$Mouse.ID)]
    X0  <- as.data.frame(covar_mat[ids_common, , drop = FALSE])
    X0  <- X0[, setdiff(colnames(X0), trait), drop = FALSE]
    
    keep <- stats::complete.cases(y0, X0)
    y0 <- y0[keep]; X0 <- X0[keep, , drop = FALSE]; Kc <- Kc[keep, keep, drop = FALSE]
    ids_keep <- rownames(Kc)
    
    if (ncol(X0)) { num <- vapply(X0, is.numeric, logical(1)); if (any(num)) X0[, num] <- lapply(X0[, num, drop=FALSE], function(x) as.numeric(scale(x, TRUE, TRUE))) }
    
    # Null model
    df0 <- cbind.data.frame(ID = factor(ids_keep, levels = ids_keep), y = y0, X0)
    rhs <- if (length(setdiff(colnames(df0), c("y","ID")))) paste(setdiff(colnames(df0), c("y","ID")), collapse = " + ") else "1"
    fml0 <- stats::as.formula(paste("y ~", rhs, "+ (1|ID)"))
    ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = maxfun))
    fit0 <- lme4qtl::relmatLmer(fml0, data = df0, relmat = list(ID = Kc), control = ctrl)
    
    sg2 <- as.numeric(lme4::VarCorr(fit0)$ID)
    se2 <- sigma(fit0)^2
    V <- sg2 * Kc + se2 * Matrix::Diagonal(nrow(Kc)); V <- (V + Matrix::t(V))/2
    eigmin <- min(eigen(as.matrix(V), symmetric = TRUE, only.values = TRUE)$values)
    if (eigmin < 1e-8) diag(V) <- diag(V) + nugget
    
    L <- chol(as.matrix(V))
    Linv <- function(a) backsolve(L, a, transpose = TRUE)
    y_t <- Linv(y0)
    Xcov <- if (ncol(X0)) cbind(Intercept = 1, X0) else matrix(1, nrow = length(y_t), dimnames = list(NULL, "Intercept"))
    X_t <- apply(Xcov, 2, Linv); if (is.vector(X_t)) X_t <- matrix(X_t, ncol = 1)
    XtX_inv <- solve(crossprod(X_t)); Xty <- crossprod(X_t, y_t)
    beta_cov <- XtX_inv %*% Xty; y_res <- y_t - X_t %*% beta_cov
    resid_against_Xt <- function(v_t) v_t - X_t %*% (XtX_inv %*% crossprod(X_t, v_t))
    sd_y <- stats::sd(y0, na.rm = TRUE)
    
    snps_cc <- markers[snp_idx]; snps_exist <- snps_cc[snps_cc %in% colnames(G)]
    G_cc <- if (length(snps_exist)) G[ids_keep, snps_exist, drop = FALSE] else NULL
    if (!is.null(G_cc) && anyNA(G_cc)) {
      cm <- colMeans(G_cc, na.rm = TRUE)
      for (k in seq_along(cm)) { idx_na <- which(is.na(G_cc[, k])); if (length(idx_na)) G_cc[idx_na, k] <- cm[k] }
    }
    maf_vec <- rep(NA_real_, length(snps_cc)); names(maf_vec) <- snps_cc
    if (length(snps_exist)) { m_vec <- colMeans(G[ids_keep, snps_exist, drop = FALSE], na.rm = TRUE); maf_exist <- pmin(m_vec/2, 1 - m_vec/2); maf_vec[names(maf_exist)] <- maf_exist }
    
    for (j in snp_idx) {
      s <- markers[j]; .i <- .i + 1L; if (verbose) cat(sprintf("[%d/%d] %s (chr %s)\n", .i, N, s, cc))
      if (!(s %in% colnames(G))) { res[[j]] <- data.frame(snp=s, chr=cc, maf=NA, beta=NA, se=NA, z=NA, p=NA, beta_orig=NA, se_orig=NA, q=NA, note="snp_not_in_G"); next }
      maf <- maf_vec[s]; if (is.na(maf) || maf < maf_min) { res[[j]] <- data.frame(snp=s, chr=cc, maf=maf, beta=NA, se=NA, z=NA, p=NA, beta_orig=NA, se_orig=NA, q=NA, note="low_var_or_maf"); next }
      g0 <- as.numeric(G_cc[, s]); g_c <- scale(g0, TRUE, FALSE)[,1]
      g_t <- Linv(g_c); g_res <- resid_against_Xt(g_t)
      Sgg <- as.numeric(crossprod(g_res)); if (Sgg <= 0) { res[[j]] <- data.frame(snp=s, chr=cc, maf=maf, beta=NA, se=NA, z=NA, p=NA, beta_orig=NA, se_orig=NA, q=NA, note="degenerate_g"); next }
      b <- as.numeric(crossprod(g_res, y_res)) / Sgg
      e <- y_res - g_res * b; n <- length(y_res); pX <- ncol(X_t); df <- max(1L, n - (pX + 1L))
      s2 <- as.numeric(crossprod(e)) / df; se <- sqrt(s2 / Sgg)
      z <- b / se; p <- 2 * stats::pnorm(-abs(z))
      res[[j]] <- data.frame(snp=s, chr=cc, maf=maf, beta=b / sd_y, se=se / sd_y, z=z, p=p, beta_orig=b, se_orig=se, q=NA, note=NA_character_, stringsAsFactors = FALSE)
    }
  }
  out <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
  if (!all(is.na(out$p))) out$q <- stats::p.adjust(out$p, method = "BH") else out$q <- NA_real_
  rownames(out) <- NULL
  out
}

# ------------------------- Trait effects helpers ----------------------------
# Standardise and subset columns for a given trait from a named list
get_trait_effects <- function(res_fast_all, trait, markers = NULL) {
  x <- res_fast_all[[trait]]
  if (is.null(x)) return(data.table::data.table())
  dt <- data.table::as.data.table(x)
  if (!"snp" %in% names(dt) && "marker" %in% names(dt)) dt[, snp := marker]
  need <- c("snp","chr","maf","beta_orig","se_orig")
  miss <- setdiff(need, names(dt)); if (length(miss)) stop(sprintf("Trait '%s' missing columns: %s", trait, paste(miss, collapse = ", ")))
  if (!is.null(markers)) dt <- dt[snp %in% markers]
  dt[, .(snp, chr, maf, beta_orig, se_orig)]
}

# ----------------------------- MR utilities ---------------------------------
# Simple IVW plot from fast-loco outputs
run_mr_from_fastloco <- function(exposure_trait, outcome_trait, res_fast_all, markers = NULL, figure = TRUE) {
  exp_dt <- get_trait_effects(res_fast_all, exposure_trait, markers)[, .(snp, chr_exp = chr, maf_exp = maf, bhat_exp = beta_orig, se_exp = se_orig)]
  out_dt <- get_trait_effects(res_fast_all, outcome_trait, markers)[, .(snp, chr_out = chr, maf_out = maf, bhat_out = beta_orig, se_out = se_orig)]
  dat <- merge(exp_dt, out_dt, by = "snp"); data.table::setDT(dat); dat <- unique(dat, by = "snp")
  dat <- dat[stats::complete.cases(bhat_exp, bhat_out, se_out)]
  if (nrow(dat) == 0L) return(NULL)
  flip <- dat$bhat_exp < 0; dat[, bhat_exp := abs(bhat_exp)]; dat[flip, bhat_out := -bhat_out]
  model <- stats::lm(bhat_out ~ 0 + bhat_exp, weights = 1/dat$se_out^2, data = dat)
  sm <- summary(model); beta <- sm$coef[1,1]; se <- sm$coef[1,2]; p <- sm$coef[1,4]
  if (isTRUE(figure)) {
    print(ggplot2::ggplot(dat, ggplot2::aes(bhat_exp, bhat_out)) +
            ggplot2::geom_point(alpha = 0.35) +
            ggplot2::geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, linewidth = 0.9) +
            ggplot2::labs(x = sprintf("β (Exposure: %s)", exposure_trait), y = sprintf("β (Outcome: %s)", outcome_trait),
                          title = sprintf("MR: %s → %s (IVW)", exposure_trait, outcome_trait),
                          subtitle = sprintf("β=%.4f ± %.4f (p=%.3g) | IVs=%d", beta, se, p, nrow(dat))) +
            ggplot2::theme_minimal(base_size = 13))
  }
  list(exposure = exposure_trait, outcome = outcome_trait, n_snp = nrow(dat), beta_mr = beta, se_mr = se, p_mr = p, dat = dat, model = model)
}

# MR with Steiger filtering (rewritten to avoid undefined helpers)
run_mr_with_steiger_markers <- function(
    exp_trait, out_trait, res_fast_all, markers_iv, N = 800,
    sf_p_threshold = 0.05, sf_use_abs = TRUE, figure = FALSE
) {
  exp_dt <- get_trait_effects(res_fast_all, exp_trait, markers_iv)[, .(id = snp, chr_exp = chr, maf_exp = maf, bhat_exp = beta_orig, se_exp = se_orig)]
  out_dt <- get_trait_effects(res_fast_all, out_trait, markers_iv)[, .(id = snp, chr_out = chr, maf_out = maf, bhat_out = beta_orig, se_out = se_orig)]
  dat <- merge(exp_dt, out_dt, by = "id"); data.table::setDT(dat)
  dat <- dat[stats::complete.cases(bhat_exp, bhat_out, se_exp, se_out, maf_exp, maf_out)]
  nsnp0 <- nrow(dat); if (nsnp0 < 3L) return(data.table::data.table(exp = exp_trait, out = out_trait, nsnp_pre = nsnp0, nsnp_post = 0, n_removed = NA_integer_, frac_removed = NA_real_, beta_pre = NA_real_, se_pre = NA_real_, p_pre = NA_real_, cochran_Q_pre = NA_real_, cochran_p_pre = NA_real_, hausman_z_pre = NA_real_, hausman_p_pre = NA_real_, beta_post = NA_real_, se_post = NA_real_, p_post = NA_real_, cochran_Q_post = NA_real_, cochran_p_post = NA_real_, hausman_z_post = NA_real_, hausman_p_post = NA_real_, sf_p_threshold = sf_p_threshold, Notes = sprintf("[Info] %s→%s: <3 IVs after gating (%d).", exp_trait, out_trait, nsnp0)))
  flip <- dat$bhat_exp < 0; dat$bhat_exp[flip] <- -dat$bhat_exp[flip]; dat$bhat_out[flip] <- -dat$bhat_out[flip]
  w <- 1/dat$se_out^2
  fit_ivw <- stats::lm(bhat_out ~ 0 + bhat_exp, weights = w, data = dat)
  sm <- summary(fit_ivw); beta_ivw <- sm$coef[1,1]; se_ivw <- sm$coef[1,2]; p_ivw <- sm$coef[1,4]
  Q <- sum(w * (dat$bhat_out - beta_ivw * dat$bhat_exp)^2); df <- nrow(dat) - 1; pQ <- stats::pchisq(Q, df, lower.tail = FALSE)
  fit_egger <- stats::lm(bhat_out ~ bhat_exp, weights = w, data = dat)
  beta_e <- coef(fit_egger)[["bhat_exp"]]; se_e <- summary(fit_egger)$coef["bhat_exp","Std. Error"]
  diff_beta <- beta_ivw - beta_e; se_diff <- sqrt((1 / sum(w * dat$bhat_exp^2)) + se_e^2)
  hz <- diff_beta / se_diff; hp <- 2 * stats::pnorm(-abs(hz))
  
  # Steiger per-IV
  dat <- dat %>% dplyr::mutate(
    R2_exp = 2 * maf_exp * (1 - maf_exp) * bhat_exp^2 / (2 * maf_exp * (1 - maf_exp) * bhat_exp^2 + se_exp^2 * N),
    R2_out = 2 * maf_out * (1 - maf_out) * bhat_out^2 / (2 * maf_out * (1 - maf_out) * bhat_out^2 + se_out^2 * N)
  )
  r_exp <- sqrt(pmax(dat$R2_exp, 0)); r_out <- sqrt(pmax(dat$R2_out, 0)); if (sf_use_abs) { r_exp <- abs(r_exp); r_out <- abs(r_out) }
  se_z <- sqrt(2 / pmax(N - 3, 1)); z_sf <- (atanh(pmin(r_exp, 0.999999)) - atanh(pmin(r_out, 0.999999))) / se_z
  dat$steiger_p <- 2 * stats::pnorm(-abs(z_sf))
  keep <- with(dat, R2_exp >= R2_out & steiger_p <= sf_p_threshold)
  dat_sf <- dat[keep, ]; nsnp1 <- nrow(dat_sf); n_removed <- nsnp0 - nsnp1; frac_removed <- n_removed / nsnp0
  
  res_post <- if (nsnp1 >= 3L) {
    w1 <- 1/dat_sf$se_out^2; fit1 <- stats::lm(bhat_out ~ 0 + bhat_exp, weights = w1, data = dat_sf)
    sm1 <- summary(fit1); list(beta = sm1$coef[1,1], se = sm1$coef[1,2], p = sm1$coef[1,4],
                               Q = sum(w1 * (dat_sf$bhat_out - sm1$coef[1,1] * dat_sf$bhat_exp)^2),
                               Qp = stats::pchisq(sum(w1 * (dat_sf$bhat_out - sm1$coef[1,1] * dat_sf$bhat_exp)^2), nrow(dat_sf)-1, lower.tail = FALSE),
                               hz = hz, hp = hp)
  } else list(beta = NA_real_, se = NA_real_, p = NA_real_, Q = NA_real_, Qp = NA_real_, hz = NA_real_, hp = NA_real_)
  
  if (isTRUE(figure)) {
    dat$.kept <- keep
    p1 <- ggplot2::ggplot(dat, ggplot2::aes(bhat_exp, bhat_out)) +
      ggplot2::geom_point(ggplot2::aes(color = .kept, size = -log10(steiger_p)), alpha = 0.6) +
      ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "steelblue3"), labels = c("Dropped (SF)", "Kept (SF)"), name = NULL) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, color = "blue", linewidth = 1.1) +
      ggplot2::geom_smooth(data = dat_sf, method = "lm", formula = y ~ 0 + x, se = FALSE, color = "red", linewidth = 1.1) +
      ggplot2::labs(title = sprintf("MR + Steiger: %s → %s", exp_trait, out_trait), x = sprintf("Beta (%s)", exp_trait), y = sprintf("Beta (%s)", out_trait)) +
      ggplot2::theme_minimal(base_size = 13)
    print(p1)
  }
  
  data.table::data.table(
    exp = exp_trait, out = out_trait,
    nsnp_pre = nsnp0, nsnp_post = nsnp1,
    n_removed = n_removed, frac_removed = frac_removed,
    beta_pre = beta_ivw, se_pre = se_ivw, p_pre = p_ivw,
    cochran_Q_pre = Q, cochran_p_pre = pQ,
    hausman_z_pre = hz, hausman_p_pre = hp,
    beta_post = res_post$beta, se_post = res_post$se, p_post = res_post$p,
    cochran_Q_post = res_post$Q, cochran_p_post = res_post$Qp,
    hausman_z_post = res_post$hz, hausman_p_post = res_post$hp,
    sf_p_threshold = sf_p_threshold
  )
}

# Normalise MR betas/SEs to per-SD units using phenotype SDs
normalize_mr_cols <- function(mr_tbl, phen, beta_col, se_col, p_col = NULL, suffix = "_norm", overwrite_fdr = FALSE) {
  stopifnot(all(c("exp","out") %in% names(mr_tbl)))
  if (!(beta_col %in% names(mr_tbl))) stop(sprintf("Column '%s' not found.", beta_col))
  if (!(se_col   %in% names(mr_tbl))) stop(sprintf("Column '%s' not found.", se_col))
  if (!is.null(p_col) && !(p_col %in% names(mr_tbl))) stop(sprintf("p_col '%s' not found.", p_col))
  
  phen_df <- as.data.frame(phen, stringsAsFactors = FALSE)
  num_idx <- vapply(phen_df, is.numeric, logical(1))
  if (!any(num_idx)) stop("`phen` has no numeric columns.")
  sd_lookup <- data.frame(trait = names(phen_df)[num_idx], sd = sapply(phen_df[num_idx], function(x) stats::sd(x, na.rm = TRUE)))
  sd_exp <- setNames(sd_lookup$sd, sd_lookup$trait); sd_out <- sd_exp
  
  scale_c <- mapply(function(e, o) {
    sde <- sd_exp[[e]]; sdo <- sd_out[[o]]; if (isTRUE(is.finite(sde) && is.finite(sdo) && sdo > 0)) sde/sdo else NA_real_
  }, mr_tbl$exp, mr_tbl$out)
  
  beta_norm_name <- paste0(beta_col, suffix); se_norm_name <- paste0(se_col, suffix)
  mr_tbl[[beta_norm_name]] <- mr_tbl[[beta_col]] * scale_c
  mr_tbl[[se_norm_name]]   <- mr_tbl[[se_col]]   * scale_c
  
  if (!is.null(p_col)) {
    fdr_target <- sub("^p", "fdr", p_col)
    fdr_name <- if (overwrite_fdr || !(fdr_target %in% names(mr_tbl))) fdr_target else paste0(fdr_target, "_recalc")
    mr_tbl[[fdr_name]] <- ifelse(is.na(mr_tbl[[p_col]]), NA_real_, stats::p.adjust(mr_tbl[[p_col]], method = "fdr"))
  }
  mr_tbl
}

# ------------------- Mouse → Human (OpenGWAS) MR API ------------------------
# OpenGWAS helpers
assert_gwas_ok <- function(id) {
  gi <- tryCatch(ieugwasr::gwasinfo(id), error = function(e) e)
  if (inherits(gi, "error") || is.null(gi) || !nrow(gi)) stop(sprintf("OpenGWAS id '%s' not accessible (check ID/JWT).", id))
  invisible(TRUE)
}

safe_extract_instruments <- function(id, pval_thresh, r2, kb, tries = 2) {
  for (k in seq_len(tries)) {
    out <- tryCatch(TwoSampleMR::extract_instruments(outcomes = id, p1 = pval_thresh, clump = TRUE, r2 = r2, kb = kb), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    if (k == tries) {
      msg <- conditionMessage(out)
      if (grepl("ieugwasr::tophits.*argument is of length zero", msg)) {
        stop(sprintf("OpenGWAS empty body for '%s' (token/rate-limit or p-value too strict).", id))
      } else stop(sprintf("Failed to extract instruments for '%s': %s", id, msg))
    }
    Sys.sleep(1)
  }
}

is_log_odds <- function(gwas_id) {
  meta <- tryCatch(ieugwasr::gwasinfo(gwas_id), error = function(e) NULL)
  if (is.null(meta)) return(FALSE)
  any(grepl("log.?odds|log\\s*OR|log-?odds", tolower(meta$unit[1])))
}

estimate_sd_from_ss <- function(dat, side = c("exposure","outcome"), gwas_id = NULL) {
  side <- match.arg(side)
  se  <- dat[[paste0("se.", side)]]
  eaf <- dat[[paste0("eaf.", side)]]
  N   <- dat[[paste0("samplesize.", side)]]
  if (all(is.na(N)) && !is.null(gwas_id)) {
    meta <- tryCatch(ieugwasr::gwasinfo(gwas_id), error = function(e) NULL)
    if (!is.null(meta) && "n" %in% names(meta) && !is.na(meta$n[1])) N <- rep(meta$n[1], nrow(dat))
  }
  ok <- is.finite(se) & is.finite(eaf) & eaf>0 & eaf<1 & is.finite(N) & N>0
  if (!any(ok)) return(NA_real_)
  stats::median(se[ok] * sqrt(N[ok]) * sqrt(2 * eaf[ok] * (1 - eaf[ok])), na.rm = TRUE)
}

human_mr_per_sd <- function(exp_id, out_id, pval_thresh = 5e-8, clump_r2 = 0.001, clump_kb = 10000, harmonise_action = 2, methods = c("mr_ivw","mr_egger_regression","mr_weighted_median")) {
  assert_gwas_ok(exp_id); assert_gwas_ok(out_id)
  exp_dat <- safe_extract_instruments(exp_id, pval_thresh, clump_r2, clump_kb); if (!nrow(exp_dat)) stop("No IVs at current thresholds.")
  out_dat <- TwoSampleMR::extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id); if (!nrow(out_dat)) stop("No overlapping SNPs in outcome for exposure IVs.")
  hdat <- TwoSampleMR::harmonise_data(exp_dat, out_dat, action = harmonise_action)
  exp_is_lo <- is_log_odds(exp_id); out_is_lo <- is_log_odds(out_id)
  sd_x <- if (!exp_is_lo) estimate_sd_from_ss(exp_dat, "exposure", exp_id) else NA_real_
  sd_y <- if (!out_is_lo) estimate_sd_from_ss(out_dat, "outcome",  out_id) else NA_real_
  if (is.finite(sd_x)) { hdat$beta.exposure <- hdat$beta.exposure/sd_x; hdat$se.exposure <- hdat$se.exposure/sd_x; hdat$exposure <- paste0(hdat$exposure," (per SD)") }
  if (is.finite(sd_y)) { hdat$beta.outcome  <- hdat$beta.outcome /sd_y; hdat$se.outcome  <- hdat$se.outcome /sd_y; hdat$outcome  <- paste0(hdat$outcome, " (per SD)") }
  mr_res <- TwoSampleMR::mr(hdat, method_list = methods)
  list(harmonised = hdat, mr_results = mr_res, sd_exposure_est = sd_x, sd_outcome_est = sd_y)
}

infer_label_from_id <- function(gwas_id, provided_label = NULL) {
  if (!is.null(provided_label) && nzchar(provided_label)) return(provided_label)
  meta <- tryCatch(ieugwasr::gwasinfo(gwas_id), error = function(e) NULL)
  if (!is.null(meta) && "trait" %in% names(meta) && nzchar(meta$trait[1])) return(meta$trait[1])
  gwas_id
}

pick_mice_effects <- function(mrow) {
  beta_col <- if ("beta_post_norm" %in% names(mrow)) "beta_post_norm" else if ("beta_norm" %in% names(mrow)) "beta_norm" else "beta"
  se_col   <- if ("se_post_norm" %in% names(mrow)) "se_post_norm" else if ("se_norm" %in% names(mrow)) "se_norm" else "se"
  list(beta = as.numeric(mrow[[beta_col]]), se = as.numeric(mrow[[se_col]]))
}

scaling_label <- function(sd_x, sd_y) {
  if (is.finite(sd_x) && is.finite(sd_y)) "IVW (per SD outcome per SD exposure)"
  else if (is.finite(sd_x))               "IVW (outcome raw units per SD exposure)"
  else if (is.finite(sd_y))               "IVW (per SD outcome per 1 unit exposure)"
  else                                    "IVW (raw units)"
}

# S3 container
new_mousehumanMR <- function(mapping, mice_row, human_result, summary) {
  structure(list(mapping = mapping, mice_row = mice_row, human_result = human_result, summary = summary), class = "mousehumanMR")
}

print.mousehumanMR <- function(x, ...) {
  s <- x$summary
  cat("Mouse → Human MR (single pair)\n")
  cat("Mouse:", s$mice_exp, "→", s$mice_out, "\n")
  cat("Human:", s$human_exp, "→", s$human_out, "\n")
  cat("IDs  :", s$exp_id, "→", s$out_id, "\n")
  cat("Method:", s$method, "\n")
  if (is.finite(s$human_beta) && is.finite(s$human_se)) {
    cat(sprintf("Human IVW: beta=%.4f, se=%.4f, p=%.3g, nsnp=%d\n", s$human_beta, s$human_se, s$human_p, s$human_nsnp))
  } else cat("Human IVW: NA\n")
  invisible(x)
}

plot.mousehumanMR <- function(x, ...) {
  s <- x$summary
  df <- data.frame(source = c("Human", "Mouse"), beta = c(s$human_beta, s$mice_beta), se = c(s$human_se, s$mice_se))
  df <- df[is.finite(df$beta) & is.finite(df$se), , drop = FALSE]
  df$lo <- df$beta - 1.96 * df$se; df$hi <- df$beta + 1.96 * df$se
  main_title <- paste0(s$mice_exp, " → ", s$mice_out, "   |   ", s$human_exp, " → ", s$human_out)
  y_lab <- paste0("Effect (", s$method, ")")
  ggplot2::ggplot(df, ggplot2::aes(x = source, y = beta)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.15) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = y_lab, title = main_title) +
    ggplot2::theme_classic(base_size = 12)
}

# Public API: Mouse→Human single-pair MR with explicit human OpenGWAS IDs
mouse_human_mr_onepair <- function(
    exp_mice, out_mice, human_exp_id, human_outcome_id, miceMR,
    human_exp_label = NULL, human_out_label = NULL,
    pval_thresh = 5e-8, clump_r2 = 0.001, clump_kb = 10000, harmonise_action = 2,
    methods = c("mr_ivw","mr_egger_regression","mr_weighted_median"), make_plot = TRUE
) {
  stopifnot(is.character(exp_mice), is.character(out_mice), is.character(human_exp_id), is.character(human_outcome_id))
  mrow <- miceMR %>% dplyr::filter(tolower(exp) == tolower(exp_mice), tolower(out) == tolower(out_mice))
  if (!nrow(mrow)) stop("No matching (exp, out) in miceMR.")
  mrow <- mrow[1, , drop = FALSE]
  
  human <- human_mr_per_sd(human_exp_id, human_outcome_id, pval_thresh, clump_r2, clump_kb, harmonise_action, methods)
  ivw <- human$mr_results %>% dplyr::filter(.data$method == "Inverse variance weighted")
  lbl_x <- infer_label_from_id(human_exp_id, human_exp_label)
  lbl_y <- infer_label_from_id(human_outcome_id, human_out_label)
  label <- scaling_label(human$sd_exposure_est, human$sd_outcome_est)
  me <- pick_mice_effects(mrow)
  
  summary_row <- tibble::tibble(
    mice_exp = mrow$exp, mice_out = mrow$out,
    human_exp = lbl_x, human_out = lbl_y,
    exp_id = human_exp_id, out_id = human_outcome_id,
    method = label,
    human_beta = if (nrow(ivw)) ivw$b else NA_real_,
    human_se   = if (nrow(ivw)) ivw$se else NA_real_,
    human_p    = if (nrow(ivw)) ivw$pval else NA_real_,
    human_nsnp = length(unique(human$harmonised$SNP)),
    human_sd_exp = human$sd_exposure_est,
    human_sd_out = human$sd_outcome_est,
    mice_beta = me$beta,
    mice_se   = me$se,
    mice_p    = mrow$p,
    mice_nsnp = mrow$nsnp
  )
  
  res <- new_mousehumanMR(
    mapping = tibble::tibble(mouse_exp = exp_mice, mouse_out = out_mice, human_exp = lbl_x, human_out = lbl_y, exp_id = human_exp_id, out_id = human_outcome_id),
    mice_row = mrow,
    human_result = human,
    summary = summary_row
  )
  if (isTRUE(make_plot)) res$plot <- plot(res)
  res
}

# ------------------------- Convenience wrappers -----------------------------
# Parallel wrapper to run Steiger+MR for one pair using provided resources
run_one_pair <- function(exp_trait, out_trait, scan_list, lod_thresh_iv, res_fast_all, N_sample, sf_p_thr = 0.05, figure = FALSE) {
  markers_iv <- get_markers_for_trait(exp_trait, scan_list, lod_thresh = lod_thresh_iv)
  n_iv <- length(markers_iv)
  if (n_iv <= 10L) {
    return(data.table::data.table(exp = exp_trait, out = out_trait, nsnp_pre = 0L, nsnp_post = 0L, n_removed = NA_integer_, frac_removed = NA_real_, beta_pre = NA_real_, se_pre = NA_real_, p_pre = NA_real_, cochran_Q_pre = NA_real_, cochran_p_pre = NA_real_, hausman_z_pre = NA_real_, hausman_p_pre = NA_real_, beta_post = NA_real_, se_post = NA_real_, p_post = NA_real_, cochran_Q_post = NA_real_, cochran_p_post = NA_real_, hausman_z_post = NA_real_, hausman_p_post = NA_real_, sf_p_threshold = sf_p_thr, Notes = sprintf("[Info] Exposure '%s': IVs ≤ 10 at LOD ≥ %g; skipped.", exp_trait, lod_thresh_iv)))
  }
  run_mr_with_steiger_markers(exp_trait, out_trait, res_fast_all, markers_iv, N = N_sample, sf_p_threshold = sf_p_thr, sf_use_abs = TRUE, figure = figure)
}

# ----------------------------------------------------------------------------

