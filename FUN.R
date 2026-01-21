# install.packages("devtools")
# devtools::install_github("variani/lme4qtl")
# ---- Libraries ----
library(dplyr)
library(stringr)
library(data.table)
library(parallel)
library(tibble)
library(ggplot2)
library(patchwork)
library(purrr)
library(MASS)
library(lme4qtl)
library(lme4)
library(lattice)
library(Matrix)
library(tidyr)
library(pbapply)
library(ggrepel)
library(viridis)
library(qtl2)
library(TwoSampleMR)
library(ieugwasr)
library(rlang)
# ---- Functions ----

get_one_region <- function(a, index, lod_thresh = 4) {
  r1 <- which(a[index:1] < lod_thresh)[1] - 2
  r2 <- which(a[index:length(a)] < lod_thresh)[1] - 2
  a[(index - r1):(index + r2)] <- 0
  return(list(a = a, index = index))
}

get_tophits <- function(a, lod_thresh = 4) {
  a[a < lod_thresh] <- 0
  peaks <- c()
  while(any(a >= lod_thresh)) {
    index <- which.max(a)
    res <- get_one_region(a, index, lod_thresh)
    peaks <- c(peaks, res$index)
    a <- res$a
  }
  return(peaks)
}

get_markers_for_trait <- function(trait, scan_list, lod_thresh = 4) {
  out <- scan_list[[trait]]
  if (is.null(out) || nrow(out) == 0) return(character())
  idx <- tryCatch(get_tophits(out, lod_thresh = lod_thresh),
                  error = function(e) integer(0))
  if (length(idx) == 0) character() else rownames(out)[idx]
}

fast_assoc <- function(y, x) {
  valid <- is.finite(y) & is.finite(x)
  y <- y[valid]
  x <- x[valid]
  n <- length(y)
  vx <- var(x)
  vy <- var(y)
  bhat <- cov(y, x) / vx
  ahat <- mean(y) - bhat * mean(x)
  rsq <- (bhat * vx)^2 / (vx * vy)
  fval <- rsq * (n - 2) / (1 - rsq)
  tval <- sqrt(fval)
  se <- abs(bhat / tval)
  pval <- pf(fval, 1, n - 2, lower.tail = FALSE)
  return(tibble(ahat = ahat, bhat = bhat, se = se, fval = fval, pval = pval, n = n))
}

get_effs <- function(g, phen, markers, trait) {
  # message("Processing trait: ", trait)
  # print(paste("Trait class is:", class(trait)))  # <- Add this>debug
  idx <- match(markers, rownames(g))
  if (any(is.na(idx))) {
    warning("Some markers not found in genotype matrix: ",
            paste(markers[is.na(idx)], collapse = ", "))
    markers <- markers[!is.na(idx)]
    idx <- idx[!is.na(idx)]
  }
  g_inst <- g[idx, , drop = FALSE]
  results <- mclapply(seq_len(nrow(g_inst)), function(i) {
    fast_assoc(phen[[trait]], g_inst[i, ]) %>%
      mutate(marker = rownames(g_inst)[i], trait = trait)
  }, mc.cores = 4)
  bind_rows(results)
}


simple_mr <- function(dat, exposure_label = "exposure_label", outcome_label = "outcome_label") {
  # Align directions so all exposure effects are positive
  flip <- dat$bhat_exp < 0
  dat$bhat_exp[flip] <- -dat$bhat_exp[flip]
  dat$bhat_out[flip] <- -dat$bhat_out[flip]
  
  # Fit inverse-variance weighted MR model (no intercept)
  model <- lm(bhat_out ~ 0 + bhat_exp, weights = 1 / dat$se_out^2, data = dat)
  coeffs <- summary(model)$coefficients
  
  # MR Plot using ggplot2
  print(
    ggplot(dat, aes(x = bhat_exp, y = bhat_out)) +
      geom_point(alpha = 0.3) +
      geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, color = "red") +
      labs(
        x = paste("Beta (Exposure:", exposure_label, ")"),
        y = paste("Beta (Outcome:", outcome_label, ")"),
        title = paste("Mendelian Randomization:", exposure_label, "→", outcome_label)
      ) +
      theme_minimal()
  )
  
  # Return result as tibble
  tibble(
    beta = coeffs[1, 1],
    se   = coeffs[1, 2],
    pval = coeffs[1, 4],
    nsnp = nrow(dat)
  )
}

# ---- inverse normal transform ----
inverse_normal_transform <- function(x, offset = 0.5, ties.method = "average") {
  # x: numeric vector (can contain NAs)
  # offset: small adjustment to avoid 0 or 1 quantiles (default 0.5)
  # ties.method: method for ranking ties ("average" is standard)
  
  # check numeric
  if (!is.numeric(x)) stop("x must be numeric.")
  
  # handle NAs
  not_na <- is.finite(x)
  r <- rank(x[not_na], ties.method = ties.method)
  n <- sum(not_na)
  
  # fractional ranks -> uniform quantiles
  u <- (r - offset) / (n - 2 * offset + 1)
  
  # inverse normal (Φ⁻¹)
  z <- qnorm(u)
  
  # return vector of same length
  out <- rep(NA_real_, length(x))
  out[not_na] <- z
  return(out)
}

get_lod_at_pos <- function(scan_obj, map_list, chr, pos, tol = 1e-4) {
  # check chromosome
  if (!chr %in% names(map_list)) {
    message(sprintf("[Warning] Chromosome '%s' not found in map.", chr))
    return(NA_real_)
  }
  
  map_chr <- map_list[[chr]]
  marker_name <- names(map_chr)[abs(map_chr - pos) < tol]
  
  if (length(marker_name) == 0) {
    message(sprintf("[Note] No marker found near chr %s at %.6f Mb (tol = %g).",
                    chr, pos, tol))
    return(NA_real_)
  }
  
  lod_value <- scan_obj[marker_name, 1]
  cat(sprintf("Chr %s pos %.6f Mb → marker %s, LOD = %.4f\n",
              chr, pos, marker_name, lod_value))
  return(lod_value)
}


run_loco_lmm_for_markers <- function(trait, markers, G, phen, covar_mat, K_loco, snp_chr,
                                     maf_min = 0.01, maxfun = 2e5, nugget = 1e-6,
                                     .print = TRUE) {
  N <- length(markers)
  res_list <- vector("list", N)

  for (i in seq_len(N)) {
    s <- markers[i]
    chr <- as.character(snp_chr[s])
    if (.print) cat(sprintf("Running [%d/%d]: %s (chr %s)\n", i, N, s, ifelse(is.na(chr) || !nzchar(chr), "NA", chr)))

    # default NA row helper
    na_row <- function(note_msg, maf = NA_real_) {
      data.frame(snp = s, chr = ifelse(is.na(chr) || !nzchar(chr), NA, chr),
                 maf = maf,
                 beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
                 beta_orig = NA_real_, se_orig = NA_real_,
                 note = note_msg, stringsAsFactors = FALSE)
    }

    # checks: SNP present?
    if (!(s %in% colnames(G))) {
      res_list[[i]] <- na_row("snp_not_in_G")
      next
    }
    # chr mapping?
    if (is.na(chr) || !nzchar(chr)) {
      res_list[[i]] <- na_row("no_chr_mapping")
      next
    }
    # LOCO K for chr?
    Kc <- K_loco[[chr]]
    if (is.null(Kc)) {
      res_list[[i]] <- na_row("no_K_for_chr")
      next
    }

    # Align to K order
    ids <- rownames(Kc)
    y0     <- phen[[trait]][match(ids, phen$Mouse.ID)]
    covar0 <- as.data.frame(covar_mat[ids, , drop = FALSE])

    # Drop the trait itself from covariates if present
    covar0 <- covar0[, setdiff(colnames(covar0), trait), drop = FALSE]

    # Drop 0-variance covariates
    if (ncol(covar0)) {
      keep_covar <- vapply(covar0, function(x) {
        if (is.factor(x)) nlevels(x) > 1 else var(as.numeric(x), na.rm = TRUE) > 0
      }, logical(1))
      covar0 <- covar0[, keep_covar, drop = FALSE]
    }

    # Keep rows with complete y & covariates
    keep_row <- complete.cases(y0, covar0)
    y0     <- y0[keep_row]
    covar0 <- covar0[keep_row, , drop = FALSE]
    Kc     <- Kc[keep_row, keep_row, drop = FALSE]
    ids_keep <- ids[keep_row]

    # Pull genotype and align, then per-SNP NA filter
    g0 <- as.numeric(G[ids_keep, s])
    if (all(is.na(g0))) {
      res_list[[i]] <- na_row("g_all_NA")
      next
    }
    keep_snp <- complete.cases(g0)
    if (!all(keep_snp)) {
      g0      <- g0[keep_snp]
      y0      <- y0[keep_snp]
      covar0  <- covar0[keep_snp, , drop = FALSE]
      Kc      <- Kc[keep_snp, keep_snp, drop = FALSE]
      ids_keep <- ids_keep[keep_snp]
    }

    # MAF / variance check
    m   <- mean(g0, na.rm = TRUE)                 # dosage mean in [0,2]
    maf <- min(m/2, 1 - m/2)
    if (is.na(maf) || maf < maf_min || sd(g0) == 0) {
      res_list[[i]] <- na_row("low_var_or_maf", maf = maf)
      next
    }

    # Scale y for stability (we’ll back-transform), scale numeric covariates; center g
    y_scaled <- scale(y0, center = TRUE, scale = TRUE)[,1]
    sd_y_orig <- sd(y0, na.rm = TRUE)

    if (ncol(covar0)) {
      num_cols <- vapply(covar0, is.numeric, logical(1))
      if (any(num_cols)) {
        covar0[, num_cols] <- lapply(covar0[, num_cols, drop = FALSE],
                                     function(x) as.numeric(scale(x, TRUE, TRUE)))
      }
    }
    g_center <- scale(g0, center = TRUE, scale = FALSE)[,1]

    # Build df, drop collapsed factors
    df <- cbind.data.frame(ID = factor(ids_keep, levels = rownames(Kc)),
                           y = y_scaled, covar0, g = g_center, row.names = NULL)
    drop_1lvl <- vapply(df, function(x) is.factor(x) && nlevels(x) < 2, logical(1))
    if (any(drop_1lvl)) df <- df[, !drop_1lvl, drop = FALSE]

    # Positive-definite nudge if necessary
    eigmin <- tryCatch(min(eigen((Kc + t(Kc))/2, symmetric = TRUE, only.values = TRUE)$values),
                       error = function(e) NA_real_)
    if (!is.na(eigmin) && eigmin < 1e-8) diag(Kc) <- diag(Kc) + nugget

    # Build formula: y ~ g + covars + (1|ID)
    cov_terms <- setdiff(colnames(df), c("y","ID","g"))
    rhs <- if (length(cov_terms)) paste(cov_terms, collapse = " + ") else "1"
    fml <- as.formula(paste("y ~ g +", rhs, "+ (1|ID)"))
    ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = maxfun))

    fit <- tryCatch(
      relmatLmer(fml, data = df, relmat = list(ID = Kc), control = ctrl),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      res_list[[i]] <- na_row(paste0("fit_error: ", fit$message), maf = maf)
      next
    }

    # Extract effects on SCALED y
    be <- tryCatch(lme4::fixef(fit)["g"], error = function(e) NA_real_)
    se <- tryCatch(sqrt(diag(vcov(fit)))["g"], error = function(e) NA_real_)
    z  <- if (is.finite(be) && is.finite(se) && se > 0) as.numeric(be/se) else NA_real_
    p  <- if (is.finite(z)) 2*pnorm(-abs(z)) else NA_real_

    # Back-transform to original phenotype scale
    beta_orig <- if (is.finite(be)) be * sd_y_orig else NA_real_
    se_orig   <- if (is.finite(se)) se * sd_y_orig else NA_real_

    res_list[[i]] <- data.frame(
      snp = s, chr = chr, maf = maf,
      beta = be, se = se, z = z, p = p,
      beta_orig = beta_orig, se_orig = se_orig,
      note = NA_character_, stringsAsFactors = FALSE
    )
  }

  ans <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  if (!is.null(ans$p)) ans$q <- p.adjust(ans$p, method = "BH")
  # already in the exact same order as `markers` by construction
  ans
}

run_fast_loco_lmm <- function(trait, markers, G, phen, covar_mat, K_loco, snp_chr,
                              maf_min = 0.01, maxfun = 2e5, nugget = 1e-6,
                              .print = TRUE) {
  
  N <- length(markers)
  res_list <- vector("list", N)
  
  # index: chr -> positions of markers in the input order
  chr_vec    <- as.character(snp_chr[markers])
  idx_by_chr <- split(seq_along(markers), chr_vec)
  
  .i <- 0L
  
  for (cc in names(idx_by_chr)) {
    if (is.na(cc) || !nzchar(cc)) {
      for (j in idx_by_chr[[cc]]) {
        s <- markers[j]
        if (.print) cat(sprintf("Running [%d/%d]: %s (chr NA)\n", { .i <- .i+1L; .i }, N, s))
        res_list[[j]] <- data.frame(snp=s, chr=NA, maf=NA,
                                    beta=NA, se=NA, z=NA, p=NA,
                                    beta_orig=NA, se_orig=NA,
                                    note="no_chr_mapping", stringsAsFactors=FALSE)
      }
      next
    }
    
    if (.print) cat("\n--- Chromosome", cc, "---\n")
    
    Kc <- K_loco[[cc]]
    if (is.null(Kc)) {
      for (j in idx_by_chr[[cc]]) {
        s <- markers[j]
        if (.print) cat(sprintf("Running [%d/%d]: %s (chr %s)\n", { .i <- .i+1L; .i }, N, s, cc))
        res_list[[j]] <- data.frame(snp=s, chr=cc, maf=NA,
                                    beta=NA, se=NA, z=NA, p=NA,
                                    beta_orig=NA, se_orig=NA,
                                    note="no_K_for_chr", stringsAsFactors=FALSE)
      }
      next
    }
    
    # ----- Align IDs across phen, covar, kinship -----
    idsK <- rownames(Kc)
    ids_common <- Reduce(intersect, list(idsK, phen$Mouse.ID, rownames(covar_mat)))
    Kc <- Kc[ids_common, ids_common, drop=FALSE]
    y0 <- phen[[trait]][match(ids_common, phen$Mouse.ID)]
    X0 <- as.data.frame(covar_mat[ids_common, , drop=FALSE])
    X0 <- X0[, setdiff(colnames(X0), trait), drop=FALSE]
    
    keep_row <- complete.cases(y0, X0)
    y0 <- y0[keep_row]
    X0 <- X0[keep_row, , drop=FALSE]
    Kc <- Kc[keep_row, keep_row, drop=FALSE]
    ids_keep <- rownames(Kc)
    
    if (ncol(X0)) {
      num_cols <- vapply(X0, is.numeric, logical(1))
      if (any(num_cols)) {
        X0[, num_cols] <- lapply(X0[, num_cols, drop=FALSE],
                                 function(x) as.numeric(scale(x, TRUE, TRUE)))
      }
    }
    
    # --- Null model on **unscaled y** ---
    df0  <- cbind.data.frame(ID=factor(ids_keep, levels=ids_keep), y=y0, X0)
    cov_terms <- setdiff(colnames(df0), c("y","ID"))
    rhs  <- if (length(cov_terms)) paste(cov_terms, collapse=" + ") else "1"
    fml0 <- as.formula(paste("y ~", rhs, "+ (1|ID)"))
    ctrl <- lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=maxfun))
    fit0 <- relmatLmer(fml0, data=df0, relmat=list(ID=Kc), control=ctrl)
    
    sg2 <- as.numeric(lme4::VarCorr(fit0)$ID)
    se2 <- sigma(fit0)^2
    
    # --- Build V, chol, whitening ---
    V <- sg2 * Kc + se2 * Diagonal(nrow(Kc))
    V <- (V + t(V))/2
    eigmin <- min(eigen(as.matrix(V), symmetric=TRUE, only.values=TRUE)$values)
    if (eigmin < 1e-8) diag(V) <- diag(V) + nugget
    
    L <- chol(as.matrix(V))
    Linv <- function(a) backsolve(L, a, transpose=TRUE)
    
    y_t <- Linv(y0)
    Xcov <- if (ncol(X0)) cbind(Intercept=1, X0) else matrix(1, nrow=length(y_t), ncol=1,
                                                             dimnames=list(NULL, "Intercept"))
    X_t <- apply(Xcov, 2, Linv)
    if (is.vector(X_t)) X_t <- matrix(X_t, ncol=1)
    
    XtX_inv <- solve(crossprod(X_t))
    Xty     <- crossprod(X_t, y_t)
    beta_cov <- XtX_inv %*% Xty
    y_res   <- y_t - X_t %*% beta_cov
    
    resid_against_Xt <- function(v_t) v_t - X_t %*% (XtX_inv %*% crossprod(X_t, v_t))
    sd_y_orig <- sd(y0, na.rm=TRUE)
    
    # --- Genotypes for this chr (aligned, imputed) ---
    snp_idx    <- idx_by_chr[[cc]]
    snps_cc    <- markers[snp_idx]
    snps_exist <- snps_cc[snps_cc %in% colnames(G)]
    G_cc <- if (length(snps_exist)) G[ids_keep, snps_exist, drop=FALSE] else NULL
    
    if (!is.null(G_cc) && anyNA(G_cc)) {
      cm <- colMeans(G_cc, na.rm=TRUE)
      for (k in seq_along(cm)) {
        idx_na <- which(is.na(G_cc[, k]))
        if (length(idx_na)) G_cc[idx_na, k] <- cm[k]
      }
    }
    
    maf_vec <- rep(NA_real_, length(snps_cc)); names(maf_vec) <- snps_cc
    if (length(snps_exist)) {
      m_vec <- colMeans(G[ids_keep, snps_exist, drop=FALSE], na.rm=TRUE)
      maf_exist <- pmin(m_vec/2, 1 - m_vec/2)
      maf_vec[names(maf_exist)] <- maf_exist
    }
    
    for (j in snp_idx) {
      s <- markers[j]
      .i <- .i + 1L
      if (.print) cat(sprintf("Running [%d/%d]: %s (chr %s)\n", .i, N, s, cc))
      
      if (!(s %in% colnames(G))) {
        res_list[[j]] <- data.frame(snp=s, chr=cc, maf=NA,
                                    beta=NA, se=NA, z=NA, p=NA,
                                    beta_orig=NA, se_orig=NA,
                                    note="snp_not_in_G", stringsAsFactors=FALSE)
        next
      }
      
      maf <- maf_vec[s]
      if (is.na(maf) || maf < maf_min) {
        res_list[[j]] <- data.frame(snp=s, chr=cc, maf=maf,
                                    beta=NA, se=NA, z=NA, p=NA,
                                    beta_orig=NA, se_orig=NA,
                                    note="low_var_or_maf", stringsAsFactors=FALSE)
        next
      }
      
      g0 <- as.numeric(G_cc[, s])
      g_c <- scale(g0, center=TRUE, scale=FALSE)[,1]
      g_t <- Linv(g_c)
      g_res <- resid_against_Xt(g_t)
      
      Sgg <- as.numeric(crossprod(g_res))
      if (Sgg <= 0) {
        res_list[[j]] <- data.frame(snp=s, chr=cc, maf=maf,
                                    beta=NA, se=NA, z=NA, p=NA,
                                    beta_orig=NA, se_orig=NA,
                                    note="degenerate_g_after_residualization", stringsAsFactors=FALSE)
        next
      }
      b  <- as.numeric(crossprod(g_res, y_res)) / Sgg
      e  <- y_res - g_res * b
      n  <- length(y_res)
      pX <- ncol(X_t)
      df <- max(1L, n - (pX + 1L))
      s2 <- as.numeric(crossprod(e)) / df
      se <- sqrt(s2 / Sgg)
      z  <- b / se
      p  <- 2 * pnorm(-abs(z))
      
      # correct: b is already on original scale; compute both versions explicitly
      beta_orig <- b
      se_orig   <- se
      beta_scaled <- b / sd_y_orig
      se_scaled   <- se / sd_y_orig
      
      res_list[[j]] <- data.frame(snp = s, chr = cc, maf = maf,
                                  beta = beta_scaled, se = se_scaled,  # scaled-y (matches slow fn)
                                  z = z, p = p,
                                  beta_orig = beta_orig, se_orig = se_orig,  # original scale
                                  note = NA_character_, stringsAsFactors = FALSE)
    }
  }
  
  # ---------- assemble & return (this was missing) ----------
  out <- rbindlist(res_list, use.names=TRUE, fill=TRUE)
  
  # if everything failed, return an empty-but-typed data.frame
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame(snp=character(), chr=character(), maf=double(),
                      beta=double(), se=double(), z=double(), p=double(),
                      beta_orig=double(), se_orig=double(),
                      q=double(), note=character(), stringsAsFactors=FALSE))
  }
  
  # BH FDR
  if (!all(is.na(out$p))) out$q <- p.adjust(out$p, method="BH") else out$q <- NA_real_
  
  # ensure same order/length as input 'markers'
  out <- out[match(markers, out$snp), ]
  rownames(out) <- NULL
  out
}


# Collect (marker, trait) pairs of top hits across all traits
collect_trait_markers <- function(scan_list, lod_thresh = 4, lod_col = 1) {
  stopifnot(is.list(scan_list))
  
  hits <- lapply(names(scan_list), function(trait) {
    out <- scan_list[[trait]]
    if (is.null(out)) return(NULL)
    if (is.vector(out)) out <- as.matrix(out)
    if (nrow(out) == 0) return(NULL)
    
    # pick LOD column (by index or name)
    if (is.character(lod_col)) {
      if (!lod_col %in% colnames(out)) return(NULL)
      lod <- out[, lod_col]
    } else {
      if (lod_col < 1 || lod_col > ncol(out)) lod_col <- 1
      lod <- out[, lod_col]
    }
    
    lod <- as.numeric(lod)
    names(lod) <- rownames(out)
    
    idx <- tryCatch(get_tophits(lod, lod_thresh = lod_thresh),
                    error = function(e) integer(0))
    if (length(idx) == 0) return(NULL)
    
    markers <- rownames(out)[idx]
    data.frame(marker = markers, trait = trait, stringsAsFactors = FALSE)
  })
  
  res <- do.call(rbind, hits)
  if (is.null(res)) {
    message("[Note] No markers passed the threshold across traits.")
    res <- data.frame(marker = character(), trait = character(), stringsAsFactors = FALSE)
  }
  return(res)
}


get_trait_effects <- function(res_fast_all, trait, markers = NULL) {
  x <- res_fast_all[[trait]]
  if (is.null(x)) {
    message(sprintf("[Note] No results for trait '%s'.", trait))
    return(data.table())
  }
  dt <- as.data.table(x)
  
  # standardize id column to 'snp' (some results may call it 'marker')
  if (!"snp" %in% names(dt) && "marker" %in% names(dt)) {
    dt[, snp := marker]
  }
  
  need <- c("snp","chr","maf","beta_orig","se_orig")
  missing_cols <- setdiff(need, names(dt))
  if (length(missing_cols)) {
    stop(sprintf("Trait '%s' results missing columns: %s",
                 trait, paste(missing_cols, collapse=", ")))
  }
  
  # optional IV filtering by markers
  if (!is.null(markers)) {
    dt <- dt[snp %in% markers]
    if (nrow(dt) == 0L) {
      message(sprintf("[Note] Trait '%s': no overlap with provided marker list.", trait))
    }
  }
  
  dt[, .(snp, chr, maf, beta_orig, se_orig)]
}

run_mr_from_fastloco <- function(
    exposure_trait,
    outcome_trait,
    res_fast_all,
    markers = NULL,     # <- OPTIONAL: restrict to this IV set
    figure  = TRUE
) {
  exp_dt <- get_trait_effects(res_fast_all, exposure_trait, markers = markers)[
    , .(snp, chr_exp = chr, maf_exp = maf,
        bhat_exp = beta_orig, se_exp = se_orig)
  ]
  
  out_dt <- get_trait_effects(res_fast_all, outcome_trait, markers = markers)[
    , .(snp, chr_out = chr, maf_out = maf,
        bhat_out = beta_orig, se_out = se_orig)
  ]
  
  # inner-join on SNP ids (your markers)
  dat <- merge(exp_dt, out_dt, by = "snp", suffixes = c("_exp","_out"))
  setDT(dat)
  dat <- unique(dat, by = "snp")
  dat <- dat[!is.na(bhat_exp) & !is.na(bhat_out) & !is.na(se_out)]
  
  if (nrow(dat) == 0L) {
    message("[Note] No overlapping IVs after filtering/merging.")
    return(NULL)
  }
  
  # align directions: force exposure betas > 0
  flip <- dat$bhat_exp < 0
  dat[,  bhat_exp := abs(bhat_exp)]
  dat[flip, bhat_out := -bhat_out]
  
  # IVW MR (no intercept), weight by outcome SE^2
  model <- lm(bhat_out ~ 0 + bhat_exp, weights = 1 / dat$se_out^2, data = dat)
  coefs <- summary(model)$coefficients
  beta_mr <- coefs[1,1]; se_mr <- coefs[1,2]; p_mr <- coefs[1,4]
  
  if (isTRUE(figure)) {
    print(
      ggplot(dat, aes(x = bhat_exp, y = bhat_out)) +
        geom_point(alpha = 0.35) +
        geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, linewidth = 0.9) +
        labs(
          x = sprintf("β (Exposure: %s)", exposure_trait),
          y = sprintf("β (Outcome: %s)", outcome_trait),
          title = sprintf("MR: %s → %s (IVW)", exposure_trait, outcome_trait),
          subtitle = sprintf("β_MR = %.4f ± %.4f  (p = %.3g)  |  IVs = %d",
                             beta_mr, se_mr, p_mr, nrow(dat))
        ) +
        theme_minimal(base_size = 13)
    )
  }
  
  list(
    exposure = exposure_trait,
    outcome  = outcome_trait,
    n_snp    = nrow(dat),
    beta_mr  = beta_mr,
    se_mr    = se_mr,
    p_mr     = p_mr,
    dat      = dat,
    model    = model
  )
}

run_mr_with_steiger_markers <- function(
    exp_trait, out_trait,
    res_fast_lmm,            # list: trait -> data.frame like your header
    markers_iv,              # IVs from get_markers_for_trait(...)
    N = 800,                 # sample size for Steiger p-values
    sf_p_threshold = 0.05,   # Steiger p-value cutoff
    sf_use_abs = TRUE,       # use |r| in Fisher z
    .figure = FALSE
) {
  # Extract & standardize columns, then gate by IVs (ids)
  exp_raw <- as.data.table(res_fast_lmm[[exp_trait]])
  out_raw <- as.data.table(res_fast_lmm[[out_trait]])
  exp_tbl <- .pick_effect_cols(exp_raw)[id %in% markers_iv]
  out_tbl <- .pick_effect_cols(out_raw)[id %in% markers_iv]
  
  # Merge exposure/outcome by id
  dat <- merge(
    exp_tbl[, .(id, chr_exp = chr, maf_exp = maf, bhat_exp = bhat, se_exp = se)],
    out_tbl[, .(id, chr_out = chr, maf_out = maf, bhat_out = bhat, se_out = se)],
    by = "id"
  )
  dat <- dat[complete.cases(dat[, .(bhat_exp, bhat_out, se_exp, se_out, maf_exp, maf_out)]), ]
  nsnp0 <- nrow(dat)
  if (nsnp0 < 3L) {
    return(data.table(
      exp = exp_trait, out = out_trait,
      nsnp_pre = nsnp0, nsnp_post = 0,
      n_removed = NA_integer_, frac_removed = NA_real_,
      beta_pre = NA_real_, se_pre = NA_real_, p_pre = NA_real_,
      cochran_Q_pre = NA_real_, cochran_p_pre = NA_real_,
      hausman_z_pre = NA_real_, hausman_p_pre = NA_real_,
      beta_post = NA_real_, se_post = NA_real_, p_post = NA_real_,
      cochran_Q_post = NA_real_, cochran_p_post = NA_real_,
      hausman_z_post = NA_real_, hausman_p_post = NA_real_,
      sf_p_threshold = sf_p_threshold,
      Notes = sprintf("[Info] %s→%s: <3 IVs after gating (%d).", exp_trait, out_trait, nsnp0)
    ))
  }
  
  # Make exposure effects positive
  flip <- dat$bhat_exp < 0
  dat$bhat_exp[flip] <- -dat$bhat_exp[flip]
  dat$bhat_out[flip] <- -dat$bhat_out[flip]
  
  # IVW + diagnostics
  compute_mr_diag <- function(d) {
    w <- 1 / d$se_out^2
    fit_ivw <- lm(bhat_out ~ 0 + bhat_exp, weights = w, data = d)
    sm <- summary(fit_ivw)
    beta_ivw <- sm$coef["bhat_exp","Estimate"]
    se_ivw   <- sm$coef["bhat_exp","Std. Error"]
    p_ivw    <- sm$coef["bhat_exp","Pr(>|t|)"]
    Q  <- sum(w * (d$bhat_out - beta_ivw * d$bhat_exp)^2)
    df <- nrow(d) - 1
    pQ <- pchisq(Q, df, lower.tail = FALSE)
    fit_egger  <- lm(bhat_out ~ bhat_exp, weights = w, data = d)
    beta_egger <- coef(fit_egger)[["bhat_exp"]]
    se_egger   <- summary(fit_egger)$coef["bhat_exp","Std. Error"]
    diff_beta <- beta_ivw - beta_egger
    se_diff   <- sqrt((1 / sum(w * d$bhat_exp^2)) + se_egger^2)
    hz        <- diff_beta / se_diff
    hp        <- 2 * pnorm(-abs(hz))
    list(beta = beta_ivw, se = se_ivw, p = p_ivw,
         cochran_Q = Q, cochran_p = pQ,
         hausman_z = hz, hausman_p = hp)
  }
  
  res_pre <- compute_mr_diag(dat)
  
  # Steiger R² and p-value per IV
  dat <- dat %>%
    mutate(
      R2_exp = 2 * maf_exp * (1 - maf_exp) * bhat_exp^2 /
        (2 * maf_exp * (1 - maf_exp) * bhat_exp^2 + se_exp^2 * N),
      R2_out = 2 * maf_out * (1 - maf_out) * bhat_out^2 /
        (2 * maf_out * (1 - maf_out) * bhat_out^2 + se_out^2 * N)
    )
  r_exp <- sqrt(pmax(dat$R2_exp, 0))
  r_out <- sqrt(pmax(dat$R2_out, 0))
  if (isTRUE(sf_use_abs)) { r_exp <- abs(r_exp); r_out <- abs(r_out) }
  se_z <- sqrt(2 / pmax(N - 3, 1))
  z_sf <- (atanh(pmin(r_exp, 0.999999)) - atanh(pmin(r_out, 0.999999))) / se_z
  p_sf <- 2 * pnorm(-abs(z_sf))
  dat$steiger_p <- p_sf
  
  # Apply Steiger filter (direction + p-value)
  keep <- with(dat, R2_exp >= R2_out & steiger_p <= sf_p_threshold)
  dat_sf <- dat[keep, ]
  nsnp1 <- nrow(dat_sf)
  n_removed <- nsnp0 - nsnp1
  frac_removed <- n_removed / nsnp0
  
  res_post <- if (nsnp1 >= 3L) compute_mr_diag(dat_sf) else
    list(beta = NA_real_, se = NA_real_, p = NA_real_,
         cochran_Q = NA_real_, cochran_p = NA_real_,
         hausman_z = NA_real_, hausman_p = NA_real_)
  
  if (.figure) {
    # Create annotation text
    lbl_stats <- sprintf(
      "Pre:  β=%.3f ± %.3f (p=%.2g)\n      Q=%.2f (p=%.2g), Hausman z=%.2f (p=%.2g)\n\nPost: β=%.3f ± %.3f (p=%.2g)\n      Q=%.2f (p=%.2g), Hausman z=%.2f (p=%.2g)\nKept: %d / %d (%.0f%%)",
      res_pre$beta,  res_pre$se,  res_pre$p,
      res_pre$cochran_Q, res_pre$cochran_p, res_pre$hausman_z, res_pre$hausman_p,
      res_post$beta, res_post$se, res_post$p,
      res_post$cochran_Q, res_post$cochran_p, res_post$hausman_z, res_post$hausman_p,
      nsnp1, nsnp0, 100*nsnp1/nsnp0
    )
    
    dat$.kept <- keep
    
    xr <- range(dat$bhat_exp, na.rm = TRUE)
    yr <- range(dat$bhat_out, na.rm = TRUE)
    x_annot <- xr[1] + 0.60 * diff(xr)
    y_annot <- yr[1] + 0.95 * diff(yr)
    
    # Separate kept and dropped IVs for distinct fits
    dat_pre <- dat
    dat_post <- dat[dat$.kept == TRUE]
    
    p1 <- ggplot(dat, aes(bhat_exp, bhat_out)) +
      geom_point(aes(color = .kept, size = -log10(steiger_p)), alpha = 0.6) +
      scale_color_manual(
        values = c("FALSE" = "grey70", "TRUE" = "steelblue3"),
        labels = c("Dropped (SF)", "Kept (SF)"),
        name = NULL
      ) +
      # Add regression lines: blue (pre), red (post)
      geom_smooth(data = dat_pre, method = "lm", formula = y ~ 0 + x,
                  se = FALSE, color = "blue", size = 1.1, linetype = "solid") +
      geom_smooth(data = dat_post, method = "lm", formula = y ~ 0 + x,
                  se = FALSE, color = "red", size = 1.1, linetype = "solid") +
      labs(
        title = sprintf("MR + Steiger Filtering: %s → %s", exp_trait, out_trait),
        subtitle = "Blue: pre-Steiger (all IVs); Red: post-Steiger (kept IVs)",
        x = sprintf("Beta (%s)", exp_trait),
        y = sprintf("Beta (%s)", out_trait)
      ) +
      annotate("label", x = x_annot, y = y_annot, hjust = 0, vjust = 1,
               label = lbl_stats, label.size = 0.25, label.padding = unit(0.3, "lines"),
               size = 3.3, alpha = 0.95) +
      guides(size = guide_legend(override.aes = list(alpha = 1), title = "-log10(SF p)")) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      )
    
    print(p1)
  }
  
  data.table(
    exp = exp_trait, out = out_trait,
    nsnp_pre = nsnp0, nsnp_post = nsnp1,
    n_removed = n_removed, frac_removed = frac_removed,
    beta_pre = res_pre$beta, se_pre = res_pre$se, p_pre = res_pre$p,
    cochran_Q_pre = res_pre$cochran_Q, cochran_p_pre = res_pre$cochran_p,
    hausman_z_pre = res_pre$hausman_z, hausman_p_pre = res_pre$hausman_p,
    beta_post = res_post$beta, se_post = res_post$se, p_post = res_post$p,
    cochran_Q_post = res_post$cochran_Q, cochran_p_post = res_post$cochran_p,
    hausman_z_post = res_post$hausman_z, hausman_p_post = res_post$hausman_p,
    sf_p_threshold = sf_p_threshold,
    Notes = ifelse(nsnp1 < 3,
                   sprintf("[Info] %s→%s: only %d IVs after Steiger (p≤%.3g).",
                           exp_trait, out_trait, nsnp1, sf_p_threshold),
                   NA_character_)
  )
}

run_one_pair <- function(exp_trait, out_trait) { # sub-function for parallel run all MR_SF
  # Get exposure IVs (marker IDs at LOD ≥ threshold)
  markers_iv <- get_markers_for_trait(exp_trait, DO_800_QTL_scan, lod_thresh = lod_thresh_iv)
  n_iv <- length(markers_iv)
  
  if (n_iv <= 10L) {
    return(data.table(
      exp = exp_trait, out = out_trait,
      nsnp_pre = 0L, nsnp_post = 0L, n_removed = NA_integer_, frac_removed = NA_real_,
      beta_pre = NA_real_, se_pre = NA_real_, p_pre = NA_real_,
      cochran_Q_pre = NA_real_, cochran_p_pre = NA_real_,
      hausman_z_pre = NA_real_, hausman_p_pre = NA_real_,
      beta_post = NA_real_, se_post = NA_real_, p_post = NA_real_,
      cochran_Q_post = NA_real_, cochran_p_post = NA_real_,
      hausman_z_post = NA_real_, hausman_p_post = NA_real_,
      sf_p_threshold = sf_p_thr,
      Notes = sprintf("[Info] Exposure '%s': IVs ≤ 10 after LOD ≥ %g; skipped.", exp_trait, lod_thresh_iv)
    ))
  }
  
  # Run MR + Steiger for this pair using your working function
  run_mr_with_steiger_markers(
    exp_trait  = exp_trait,
    out_trait  = out_trait,
    res_fast_lmm = res_fast_all,       # your list of per-trait fast LMM results
    markers_iv   = markers_iv,
    N = N_sample,
    sf_p_threshold = sf_p_thr,
    .figure = FALSE                    # set TRUE if you want plots per pair
  )
}

normalize_mr_cols <- function(
    mr_tbl,
    phen,
    beta_col,              # e.g. "beta_post"
    se_col,                # e.g. "se_post"
    p_col      = NULL,     # e.g. "p_post" (optional)
    suffix     = "_norm",  # appended to beta/se colnames
    overwrite_fdr = FALSE  # if TRUE, overwrite matching FDR (e.g., fdr_post)
) {
  # --- checks
  if (!all(c("exp", "out") %in% names(mr_tbl))) {
    stop("`mr_tbl` must have columns 'exp' and 'out'.")
  }
  if (!(beta_col %in% names(mr_tbl))) stop(sprintf("Column '%s' not found.", beta_col))
  if (!(se_col   %in% names(mr_tbl))) stop(sprintf("Column '%s' not found.", se_col))
  if (!is.null(p_col) && !(p_col %in% names(mr_tbl))) {
    stop(sprintf("p_col '%s' not found in mr_tbl.", p_col))
  }
  
  # --- SD lookup from phenotype matrix/data.frame (base R)
  phen_df <- as.data.frame(phen, stringsAsFactors = FALSE)
  num_idx <- vapply(phen_df, is.numeric, logical(1))
  if (!any(num_idx)) stop("`phen` has no numeric columns.")
  
  sd_lookup <- data.frame(
    trait = names(phen_df)[num_idx],
    sd    = sapply(phen_df[num_idx], function(x) stats::sd(x, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  # --- merge SDs into MR table
  sd_exp <- sd_lookup
  names(sd_exp) <- c("exp", "sd_exp")
  sd_out <- sd_lookup
  names(sd_out) <- c("out", "sd_out")
  
  out <- merge(mr_tbl, sd_exp, by = "exp", all.x = TRUE)
  out <- merge(out,    sd_out, by = "out", all.x = TRUE)
  
  # --- compute scale = sd_exp / sd_out and normalized columns
  out$scale_c <- with(out, ifelse(is.finite(sd_exp) & is.finite(sd_out) &
                                    !is.na(sd_exp) & !is.na(sd_out) & sd_out > 0,
                                  sd_exp / sd_out, NA_real_))
  
  beta_norm_name <- paste0(beta_col, suffix)
  se_norm_name   <- paste0(se_col,   suffix)
  
  out[[beta_norm_name]] <- out[[beta_col]] * out[["scale_c"]]
  out[[se_norm_name]]   <- out[[se_col]]   * out[["scale_c"]]
  
  # --- recompute FDR from p_col (optional)
  if (!is.null(p_col)) {
    fdr_target <- sub("^p", "fdr", p_col)
    fdr_name <- if (overwrite_fdr || !(fdr_target %in% names(out))) {
      fdr_target
    } else {
      paste0(fdr_target, "_recalc")
    }
    out[[fdr_name]] <- ifelse(is.na(out[[p_col]]), NA_real_,
                              stats::p.adjust(out[[p_col]], method = "fdr"))
  }
  
  out
}


# single-pair mouse→human MR with OpenGWAS IDs and optional labels.
# ------------------------- utilities --------------------------
assert_gwas_ok <- function(id) {
  gi <- tryCatch(ieugwasr::gwasinfo(id), error = function(e) e)
  if (inherits(gi, "error") || is.null(gi) || !nrow(gi)) {
    stop(sprintf("OpenGWAS id '%s' not accessible (check ID or JWT token).", id))
  }
  invisible(TRUE)
}

safe_extract_instruments <- function(id, pval_thresh, r2, kb, tries = 2) {
  for (k in seq_len(tries)) {
    out <- tryCatch(
      TwoSampleMR::extract_instruments(outcomes = id, p1 = pval_thresh, clump = TRUE, r2 = r2, kb = kb),
      error = function(e) e
    )
    if (!inherits(out, "error")) return(out)
    if (k == tries) {
      msg <- conditionMessage(out)
      if (grepl("ieugwasr::tophits.*argument is of length zero", msg)) {
        stop(sprintf(
          "OpenGWAS returned an empty body for '%s'. Possible token/rate-limit or stringent p-value.\n- Confirm the ID\n- Ensure JWT via ieugwasr::set_opengwas_jwt()\n- Try a looser pval_thresh.",
          id
        ))
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

human_mr_per_sd <- function(exp_id, out_id,
                            pval_thresh = 5e-8, clump_r2 = 0.001, clump_kb = 10000,
                            harmonise_action = 2,
                            methods = c("mr_ivw","mr_egger_regression","mr_weighted_median")) {
  assert_gwas_ok(exp_id); assert_gwas_ok(out_id)
  
  exp_dat <- safe_extract_instruments(exp_id, pval_thresh, clump_r2, clump_kb)
  if (!nrow(exp_dat)) stop("No IVs at current thresholds.")
  
  out_dat <- TwoSampleMR::extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id)
  if (!nrow(out_dat)) stop("No overlapping SNPs in outcome for exposure IVs.")
  
  hdat <- TwoSampleMR::harmonise_data(exp_dat, out_dat, action = harmonise_action)
  
  exp_is_lo <- is_log_odds(exp_id); out_is_lo <- is_log_odds(out_id)
  sd_x <- if (!exp_is_lo) estimate_sd_from_ss(exp_dat, "exposure", exp_id) else NA_real_
  sd_y <- if (!out_is_lo) estimate_sd_from_ss(out_dat, "outcome",  out_id) else NA_real_
  
  if (is.finite(sd_x)) {
    hdat$beta.exposure <- hdat$beta.exposure/sd_x
    hdat$se.exposure   <- hdat$se.exposure/sd_x
    hdat$exposure      <- paste0(hdat$exposure," (per SD)")
  }
  if (is.finite(sd_y)) {
    hdat$beta.outcome <- hdat$beta.outcome/sd_y
    hdat$se.outcome   <- hdat$se.outcome/sd_y
    hdat$outcome      <- paste0(hdat$outcome," (per SD)")
  }
  
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
  beta_col <- if ("beta_post_norm" %in% names(mrow)) "beta_post_norm"
  else if ("beta_norm" %in% names(mrow)) "beta_norm" else "beta"
  se_col   <- if ("se_post_norm" %in% names(mrow)) "se_post_norm"
  else if ("se_norm"   %in% names(mrow)) "se_norm" else "se"
  list(beta = as.numeric(mrow[[beta_col]]), se = as.numeric(mrow[[se_col]]))
}

scaling_label <- function(sd_x, sd_y) {
  if (is.finite(sd_x) && is.finite(sd_y)) "IVW (per SD outcome per SD exposure)"
  else if (is.finite(sd_x))                "IVW (outcome raw units per SD exposure)"
  else if (is.finite(sd_y))                "IVW (per SD outcome per 1 unit exposure)"
  else                                     "IVW (raw units)"
}

# --------------------------- S3 class --------------------------
new_mousehumanMR <- function(mapping, mice_row, human_result, summary) {
  structure(
    list(mapping = mapping, mice_row = mice_row, human_result = human_result, summary = summary),
    class = "mousehumanMR"
  )
}

print.mousehumanMR <- function(x, ...) {
  s <- x$summary
  cat("Mouse → Human MR (single pair)\n")
  cat("Mouse:", s$mice_exp, "→", s$mice_out, "\n")
  cat("Human:", s$human_exp, "→", s$human_out, "\n")
  cat("IDs  :", s$exp_id, "→", s$out_id, "\n")
  cat("Method:", s$method, "\n")
  if (is.finite(s$human_beta) && is.finite(s$human_se)) {
    z <- s$human_beta / s$human_se
    cat(sprintf("Human IVW: beta=%.4f, se=%.4f, p=%.3g, nsnp=%d\n",
                s$human_beta, s$human_se, s$human_p, s$human_nsnp))
  } else {
    cat("Human IVW: NA (check instruments/overlap)\n")
  }
  invisible(x)
}

plot.mousehumanMR <- function(x, ...) {
  s <- x$summary
  df <- data.frame(
    source = c("Human", "Mouse"),
    beta   = c(s$human_beta, s$mice_beta),
    se     = c(s$human_se,   s$mice_se)
  )
  df <- df[is.finite(df$beta) & is.finite(df$se), , drop = FALSE]
  df$lo <- df$beta - 1.96*df$se
  df$hi <- df$beta + 1.96*df$se
  
  main_title <- paste0(s$mice_exp, " → ", s$mice_out, "   |   ",
                       s$human_exp, " → ", s$human_out)
  y_lab <- paste0("Effect (", s$method, ")")
  
  ggplot2::ggplot(df, ggplot2::aes(x = source, y = beta)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.15) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = y_lab, title = main_title) +
    ggplot2::theme_classic(base_size = 12)
}

# --------------------------- main API -------------------------
#' Mouse→Human MR for a single pair (mouse traits fixed; human specified by OpenGWAS IDs)
#'
#' @param exp_mice Character. Mouse exposure trait name (matches miceMR$exp).
#' @param out_mice Character. Mouse outcome trait name (matches miceMR$out).
#' @param human_exp_id Character. OpenGWAS ID for human exposure.
#' @param human_outcome_id Character. OpenGWAS ID for human outcome.
#' @param miceMR Data.frame with columns: exp, out, p, nsnp, and mouse effects
#'        (preferably beta_post_norm, se_post_norm; fallbacks to beta_norm/se_norm or beta/se).
#' @param human_exp_label Character or NULL. Label to display for human exposure (defaults to trait from OpenGWAS).
#' @param human_out_label Character or NULL. Label to display for human outcome (defaults to trait from OpenGWAS).
#' @param pval_thresh,clump_r2,clump_kb,harmonise_action MR extraction/harmonisation settings.
#' @param methods Character vector of TwoSampleMR methods to run (default IVW/Egger/WM).
#' @param make_plot Logical. If TRUE, attaches a ggplot under $plot for convenience.
#' @return An object of class "mousehumanMR" with elements: mapping, mice_row, human_result, summary, and optionally plot.
mouse_human_mr_onepair <- function(
    exp_mice,
    out_mice,
    human_exp_id,
    human_outcome_id,
    miceMR,
    human_exp_label = NULL,
    human_out_label = NULL,
    pval_thresh = 5e-8,
    clump_r2 = 0.001,
    clump_kb = 10000,
    harmonise_action = 2,
    methods = c("mr_ivw","mr_egger_regression","mr_weighted_median"),
    make_plot = TRUE
) {
  stopifnot(is.character(exp_mice), is.character(out_mice),
            is.character(human_exp_id), is.character(human_outcome_id))
  
  # 1) locate the mouse row (first match if multiple)
  mrow <- miceMR %>%
    dplyr::filter(tolower(exp) == tolower(exp_mice),
                  tolower(out) == tolower(out_mice))
  if (!nrow(mrow)) stop("No matching (exp, out) in miceMR.")
  mrow <- mrow[1, , drop = FALSE]
  
  # 2) run human MR
  human <- human_mr_per_sd(
    exp_id = human_exp_id,
    out_id = human_outcome_id,
    pval_thresh = pval_thresh,
    clump_r2 = clump_r2,
    clump_kb = clump_kb,
    harmonise_action = harmonise_action,
    methods = methods
  )
  
  # 3) IVW row & labels
  ivw <- human$mr_results %>% dplyr::filter(.data$method == "Inverse variance weighted")
  lbl_x <- infer_label_from_id(human_exp_id, human_exp_label)
  lbl_y <- infer_label_from_id(human_outcome_id, human_out_label)
  label <- scaling_label(human$sd_exposure_est, human$sd_outcome_est)
  
  # 4) mouse effects (post_norm first)
  me <- pick_mice_effects(mrow)
  
  # 5) summary row
  summary_row <- tibble::tibble(
    mice_exp = mrow$exp,
    mice_out = mrow$out,
    human_exp = lbl_x,
    human_out = lbl_y,
    exp_id = human_exp_id,
    out_id = human_outcome_id,
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
    mapping = tibble::tibble(
      mouse_exp = exp_mice, mouse_out = out_mice,
      human_exp = lbl_x, human_out = lbl_y,
      exp_id = human_exp_id, out_id = human_outcome_id
    ),
    mice_row = mrow,
    human_result = human,
    summary = summary_row
  )
  
  if (isTRUE(make_plot)) res$plot <- plot(res)
  res
}