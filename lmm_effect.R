# NOTICE: 
# 1, the qtl2 package allele prob is based on 8 founder can not be 0/1/2
# 2, what we want is using LMM to consider covariate, Kinship and result in 
#    genotype 0/1/2 association with phenotype for particular SNPs
# 
# DO_800_QTL_scan <- scan1(
#   probs    = DO_800_allele_probs,
#   pheno    = cross_basic$pheno,
#   kinship  = DO_800_kinship_loco,
#   addcovar = cross_basic$covar
# )
# ---- FUN.R ----
source("./FUN.R")
# ---- Load Data ----
datadir <- "/Users/xd14188/Desktop/UoB/miceMR/working/data/Farber_Lab_DO_genotype"
# map <- fread(file.path(datadir, "Marker-positions/GM_markers_pmap.csv"))
load(file.path(datadir, "QTL2_RData/DO_800_QTL_scan.RData"))
# load(file.path(datadir, "QTL2_RData/DO_800_allele_probs.RData"))
load(file.path(datadir, "QTL2_RData/DO_800_kinship_loco.RData"))
load(file.path(datadir, "QTL2_RData/cross_basic_cleaned.Rdata"))

# Get the phenotype and the covariate data from the cross object
# Phenotype data
DO_800_pheno_data <- as.data.frame(cross_basic$pheno)
# Covariate data, do not use body length and body weight in Covariate
DO_800_covar_data_matrix <-
  cross_basic$covar %>%
  as.data.frame() %>%
  rownames_to_column("Mouse.ID") %>%
  transmute(
    Mouse.ID,
    sex = as.integer(sex == "M"),
    age_at_sac_days = as.numeric(age_at_sac_days),
    ngen = as.numeric(ngen)
  ) %>%
  column_to_rownames("Mouse.ID") %>%
  as.matrix(); head(DO_800_covar_data_matrix)

datadir <- "/Users/xd14188/Desktop/UoB/miceMR/working/data/Farber_Lab_DO_genotype"
geno <- fread(file.path(datadir, "Genotype/DO_800_geno_full.csv"), header = T)
phen <- fread(file.path(datadir, "Phenotype/DO_800_pheno_full_with_covars.csv"))

# ---- Preprocessing Genotype ----
markers <- rownames(DO_800_QTL_scan[[1]])
geno <- geno[match(markers, geno$Marker), ]
stopifnot(all(geno$Marker == markers))

g <- as.matrix(geno[, -1])
rownames(g) <- geno$Marker

# Align phenotype and genotype samples
common_ids <- intersect(phen$Mouse.ID, colnames(g))
phen <- phen[phen$Mouse.ID %in% common_ids, ]
g <- g[, match(phen$Mouse.ID, colnames(g))]
stopifnot(all(phen$Mouse.ID == colnames(g)))
# can not add BMI as scan results DO_800_QTL_scan do not have BMI
# phen$BMI <- phen$body_weight / (phen$body_length ^ 2) # adding BMI into phen

# Prepare the chr for the SNPs
pmap_list <- cross_basic$pmap
snp_chr <- unlist(lapply(names(pmap_list), function(chr) {
  setNames(rep(chr, length(pmap_list[[chr]])), names(pmap_list[[chr]]))
}))

#########
# Start #
#########
marker_trait_table <- collect_trait_markers(DO_800_QTL_scan, lod_thresh = 4)
head(marker_trait_table)
marker_all <- unique(marker_trait_table$marker); dim(marker_trait_table); length(marker_all)

traits <- names(DO_800_QTL_scan)  # or your subset
cores  <- 4  # max(1L, detectCores() - 1L)

res_fast_all <- mclapply(traits, function(tr) {
  message(sprintf("==> %s (%d markers)", tr, length(marker_all)))
  tryCatch(
    run_fast_loco_lmm(
      trait     = tr,
      markers   = marker_all,
      G         = t(g),
      phen      = phen,
      covar_mat = DO_800_covar_data_matrix,
      K_loco    = DO_800_kinship_loco,
      snp_chr   = snp_chr,
      .print    = FALSE    
    ),
    error = function(e) { message(sprintf("[Warn] %s failed: %s", tr, e$message)); NULL }
  )
}, mc.cores = cores)

names(res_fast_all) <- traits


head(res_fast_all)
save(res_fast_all, file = "res_lmm.rda")
# The LMM loco opt demo
# trait <- "body_weight"
# out   <- DO_800_QTL_scan[[trait]]
# th <- get_tophits(out, lod_thresh=4)
# markers <- rownames(DO_800_QTL_scan[[trait]])[th]
# effs <- get_effs(g, phen, markers, trait)
# system.time(
#   effs_lmm <- run_loco_lmm_for_markers(trait, markers, G = t(g), phen, DO_800_covar_data_matrix,
#                                        DO_800_kinship_loco, snp_chr)
# ) # 650.870s for 93 SNPs
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# investigate some particular SNP 
# from the fast_associate this SNP result in outlier result from MAF << 0.01
# we can find it and filter in both lmm and lmm_fast
# below is the sample 
# ------------------------------------------------------------------------------
# trait <- "gastroc_weight"
# idx <- 31 # which.max(abs(effs$bhat)) # 31
# effs[idx, ]
# effs_lmm[idx, ]
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# LINEAR MIXED MODEL (LMM) GWAS — "OPTIMIZATION" VERSION SUMMARY
# ------------------------------------------------------------------------------
# This function fits a linear mixed model (LMM) for each SNP using relmatLmer(),
# which *optimizes* variance components (σ_g², σ_e²) individually for every SNP.
# It is statistically rigorous but computationally slow (~seconds per SNP).
#
# MODEL:
#   y = Xβ + Zu + ε
#   u ~ N(0, σ_g² * K)        # random genetic effects based on kinship K
#   ε ~ N(0, σ_e² * I)        # residual noise
#   → Var(y) = σ_g²*K + σ_e²*I
#
# RELMATLMER ESTIMATES:
#   - Fixed effects β (covariates + SNP effect)
#   - Variance components σ_g², σ_e² via REML optimisation
#   - Uses "bobyqa" optimizer (bounded quasi-Newton) to maximize REML likelihood
#
# BASIC STEPS:
#   1. Align phenotype, genotype, and covariates to the same sample order.
#   2. Build model formula: y ~ g + covariates + (1|ID)
#   3. For each SNP:
#        • Center genotype (0/1/2 dosage)
#        • Check MAF and variance
#        • Ensure kinship Kc is positive-definite (add tiny 'nugget' if needed)
#        • Call relmatLmer() to estimate σ_g², σ_e², β, SE, z, p
#        • Back-transform β, SE to original phenotype scale
#   4. Combine all SNP results into one data.frame (same order as input)
#
# PARAMETER SETTINGS:
#   maf_min = 0.01     # skip SNPs with MAF below this threshold
#   maxfun  = 2e5      # max optimizer iterations for convergence
#   nugget  = 1e-6     # small diagonal added if K nearly singular
#   .print  = TRUE     # display progress: [i/total] SNP name and chr
#
# INTUITIVE IDEA:
#   - The LMM accounts for relatedness via kinship K, preventing false positives.
#   - relmatLmer() iteratively adjusts σ_g² and σ_e² so the model fits the data
#     best (maximum REML likelihood). This iterative "optimization" is why
#     the method is accurate but computationally heavy.
#   - Each SNP is tested in a separate LMM, so variance components are refit
#     each time, making the approach slow but statistically exact.
#
# PRACTICAL NOTES:
#   • Produces β (scaled y) and β_orig (original-scale effect size)
#   • Adds minor allele frequency (maf) to results
#   • Adds a small nugget to K to avoid singularities
#   • Outputs are in the same order and length as input 'markers'
#
# This version is the "gold standard" LMM — exact but slow.
# For large-scale scans, consider a "null model once + GLS" approach (EMMAX style).
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# WHY relmatLmer() IS SLOW AND WHY "OPTIMISATION" IS NEEDED
# ------------------------------------------------------------------------------
# Linear mixed models (LMM) include random effects to model sample relatedness:
#     y = X*beta + u + e
#   where u ~ N(0, σ_g^2 * K) and e ~ N(0, σ_e^2 * I)
#
# The unknown variance components σ_g^2 (genetic) and σ_e^2 (residual)
# cannot be solved analytically. relmatLmer() must *optimise* them by
# maximising the (restricted) likelihood over many iterations.
#
# Each optimisation step requires matrix decompositions of the n×n covariance
# matrix Σ = σ_g^2*K + σ_e^2*I  →  computationally heavy (≈ O(n^3)).
#
# Because relmatLmer() re-estimates σ_g^2 and σ_e^2 for *every SNP*,
# it performs a full optimisation loop each time → slow (~seconds per SNP).
#
# Fast GWAS tools (e.g. GEMMA, EMMAX, fastGWA) fit the null model once,
# save the estimated variance components, and then test all SNPs by
# Generalised Least Squares (GLS) without re-optimising → orders of magnitude faster.
# ------------------------------------------------------------------------------
# system.time(
# res_fast <- run_fast_loco_lmm( trait, markers, G = t(g), phen, 
#                                DO_800_covar_data_matrix, DO_800_kinship_loco, snp_chr,
#                                .print = TRUE)
# ) # 115.913s for 93 SNPs
# 
# head(res_fast)
# ------------------------------------------------------------------------------
# FAST LINEAR MIXED MODEL (LMM) GWAS — LOCO (CHROMOSOME-SPECIFIC) GLS VERSION
# ------------------------------------------------------------------------------
# This function implements an *approximate but much faster* LMM scan using
# generalized least squares (GLS) with a LOCO kinship matrix.
#
# It reuses variance components (σ_g², σ_e²) estimated from a single "null"
# model per chromosome (without SNP term) to test all SNPs on that chromosome.
# This avoids re-optimizing variance parameters for every SNP, which is the
# main time cost in full per-SNP LMM fitting.
#
# MODEL:
#   y = Xβ + gγ + u + ε
#   u ~ N(0, σ_g² * K_chr)      # random genetic effects (LOCO kinship)
#   ε ~ N(0, σ_e² * I)
#   → Var(y) = σ_g²*K_chr + σ_e²*I = V
#
# PROCEDURE:
#   1. Fit a *null* mixed model (y ~ covariates + (1|ID)) for each chromosome:
#        → Obtain σ̂_g² and σ̂_e² once per chromosome (no SNPs yet).
#   2. Compute V = σ̂_g²*K_chr + σ̂_e²*I and its Cholesky factor L = chol(V).
#   3. "Whiten" data:  y_t = L⁻¹y,  X_t = L⁻¹X
#   4. Precompute projection matrices for covariates:
#        XtX_inv = (X_t'X_t)⁻¹, etc.
#   5. For each SNP on this chromosome:
#        • mean-impute missing genotypes
#        • center g, whiten g_t = L⁻¹g
#        • residualize g_t and y_t against X_t (Frisch–Waugh–Lovell theorem)
#        • regress residualized y_t on residualized g_t
#        → analytic β̂, SE, z, p
#
# SPEEDUP IDEA:
#   - Full LMM: variance parameters re-optimized by REML for *each SNP*.
#   - Fast GLS: variance parameters fixed (from null model), only GLS regression
#     per SNP (~linear time, no iterative optimization).
#
# OUTPUT:
#   - β, SE, z, p : effect estimates on *scaled-y* scale (comparable across traits)
#   - β_orig, SE_orig : effect estimates on *original* phenotype scale
#   - maf : minor allele frequency per SNP
#   - q : BH-adjusted FDR
#   - note : annotation for low-MAF or missing SNPs
#   - results returned in same order as input 'markers'
#
# PARAMETER SETTINGS:
#   maf_min = 0.01   # minimum MAF threshold
#   nugget  = 1e-6   # stabilizer for near-singular K matrices
#   maxfun  = 2e5    # optimization limit for null model REML fit
#   .print  = TRUE   # progress output per SNP
#
# INTUITIVE IDEA:
#   - By fitting the null model once per chromosome, we treat the kinship-based
#     variance structure as fixed across SNPs in that region (a good LOCO
#     approximation).
#   - Whitening by the Cholesky of V removes sample correlations, so testing
#     each SNP reduces to a simple linear regression in whitened space.
#   - This makes the scan *tens to hundreds of times faster* than full LMM,
#     with negligible loss in accuracy for most GWAS settings.
#
# PRACTICAL NOTES:
#   • β_orig and SE_orig are on the same scale as the input phenotype.
#   • β and SE are standardized by SD(y) for comparability.
#   • Mean-imputation of genotypes is used to keep sample size consistent.
#   • Uses LOCO kinship (leave-one-chromosome-out) to avoid proximal contamination.
#
# This approach reproduces full-LMM results very closely but runs far faster.
# It’s conceptually similar to GEMMA/EMMAX/fastGWA GLS methods.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#  my next idea is to have the all tophits (from all phenotypes) to make a 
#  long vector of markers, which is are the input for all LMM fast
#  to find the beta and s.e., for different phenotype 
#  tring to run the beta calculation once 
#  as the LMM fast only slow for change chr to a new chr and 
#  very fast within the same chr, so adding the number of SNPs does not matter
#  it results in a big effect and s.e., table contain all the VI for different trait
#  THE CODE FOR RUN ALL MAKRERS CORSSING ALL TRAIT AT THE TOP OF COMMENT
# ------------------------------------------------------------------------------
