# using the simple MR to check the effects from LMM compare with 
# MR results using effects from fast_associate

exp <- "ML"
trait <- exp
out   <- DO_800_QTL_scan[[trait]]
th <- get_tophits(out, lod_thresh=4)
markers <- rownames(DO_800_QTL_scan[[trait]])[th]
effs_exp <- get_effs(g, phen, markers, trait)
system.time(
  effs_lmm_exp <- run_loco_lmm_for_markers(trait, markers, G = t(g), phen, DO_800_covar_data_matrix,
                                       DO_800_kinship_loco, snp_chr) 
) 

outcome <- "uCT_pMOI"
trait <- outcome
# out   <- DO_800_QTL_scan[[trait]]
# th <- get_tophits(out, lod_thresh=4)
# markers <- rownames(DO_800_QTL_scan[[trait]])[th]
effs_outcome <- get_effs(g, phen, markers, trait)
system.time(
  effs_lmm_outcome <- run_loco_lmm_for_markers(trait, markers, G = t(g), phen, DO_800_covar_data_matrix,
                                           DO_800_kinship_loco, snp_chr) 
) 


a <- effs_lmm_exp$beta
b <- effs_lmm_outcome$beta

plot(a, b)



## --- 1) MR with original-scale effects (recommended) ---
# Prepare exposure side
exp_dt <- as.data.table(effs_lmm_exp)[
  , .(snp, chr_exp = chr, maf_exp = maf,
      bhat_exp = beta_orig, se_exp = se_orig)
]

# Prepare outcome side
out_dt <- as.data.table(effs_lmm_outcome)[
  , .(snp, chr_out = chr, maf_out = maf,
      bhat_out = beta_orig, se_out = se_orig)
]

# Inner-join by SNP, keep one row per SNP, drop NAs needed by simple_mr()
dat <- merge(exp_dt, out_dt, by = "snp", suffixes = c("_exp_ML","_out_uCT_pMOI"))
dat <- unique(dat, by = "snp")
dat <- dat[!is.na(bhat_exp) & !is.na(bhat_out) & !is.na(se_out)]

# run MR
res_mr <- simple_mr(dat, exposure_label = "_exp_ML", outcome_label = "_out_uCT_pMOI")
print(res_mr)

