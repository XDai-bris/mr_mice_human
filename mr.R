source("./FUN.R")
load("./res_lmm.rda")
datadir <- "/Users/xd14188/Desktop/UoB/miceMR/working/data/Farber_Lab_DO_genotype"
# map <- fread(file.path(datadir, "Marker-positions/GM_markers_pmap.csv"))
load(file.path(datadir, "QTL2_RData/DO_800_QTL_scan.RData"))
# load(file.path(datadir, "QTL2_RData/DO_800_allele_probs.RData"))
load(file.path(datadir, "QTL2_RData/DO_800_kinship_loco.RData"))
load(file.path(datadir, "QTL2_RData/cross_basic_cleaned.Rdata"))

# Get the phenotype and the covariate data from the cross object
# Covariate data
DO_800_covar_data <- as.data.frame(cross_basic$covar)
# Add the "Mouse.ID" column in the covar data
DO_800_covar_data <- DO_800_covar_data %>% rownames_to_column(var="Mouse.ID")
# Remove extra columns from the covar data
DO_800_covar_data <- DO_800_covar_data[, c(1:7)]
# Phenotype data
DO_800_pheno_data <- as.data.frame(cross_basic$pheno)
# Encode sex as numeric values in the covariate data
DO_800_covar_data[,"sex"] = (DO_800_covar_data[,"sex"] == "M")
# Create a numeric matrix with the covariate data
DO_800_covar_data_matrix = apply(DO_800_covar_data,2,as.numeric)
# Make sure the rownames of the numeric matrix corresponds to mouse Ids
rownames(DO_800_covar_data_matrix) <- rownames(cross_basic$covar)
# Remove the Mouse.ID and the sac_date columns from the covariate
DO_800_covar_data_matrix <- DO_800_covar_data_matrix[,- c(1:2)]
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

# Prepare the chr for the SNPs
pmap_list <- cross_basic$pmap
snp_chr <- unlist(lapply(names(pmap_list), function(chr) {
  setNames(rep(chr, length(pmap_list[[chr]])), names(pmap_list[[chr]]))
}))

# Pick two traits from res_fast_all
exp_trait <- "uCT_pMOI"
out   <- DO_800_QTL_scan[[exp_trait]]
th <- get_tophits(out, lod_thresh=4)
markers <- rownames(DO_800_QTL_scan[[exp_trait]])[th]
out_trait <- "ML"

# Run MR (with figure)
res_mr <- run_mr_from_fastloco(
  exposure_trait = exp_trait,
  outcome_trait  = out_trait,
  res_fast_all   = res_fast_all,
  markers        = markers,   # <-- filter by these IVs
  figure         = TRUE          # turn FALSE to suppress plot
)


# View summary
str(res_mr, 1)



# run all of MR

traits <- intersect(names(res_fast_all), names(DO_800_QTL_scan))
cores  <- max(1L, detectCores() - 1L)

get_markers_for_trait <- function(trait, scan_list, lod_thresh = 4) {
  out <- scan_list[[trait]]
  if (is.null(out) || nrow(out) == 0) return(character())
  idx <- tryCatch(get_tophits(out, lod_thresh = lod_thresh),
                  error = function(e) integer(0))
  if (length(idx) == 0) character() else rownames(out)[idx]
}

res_list <- mclapply(seq_along(traits), function(i) {
  exp_trait <- traits[i]
  markers_iv <- get_markers_for_trait(exp_trait, DO_800_QTL_scan, lod_thresh = 4)
  n_iv <- length(markers_iv)
  
  # Handle missing or too few IVs (record Notes instead of skipping silently)
  if (n_iv == 0L) {
    return(data.table(
      exp   = exp_trait, out = NA_character_,
      exp_i = i, out_i = NA_integer_,
      beta  = NA_real_, se = NA_real_, p = NA_real_,
      nsnp  = NA_integer_, fdr = NA_real_,
      Notes = sprintf("[Info] Exposure '%s': no IVs (LOD ≥ 4); skipped.", exp_trait)
    ))
  }
  
  if (n_iv <= 10L) {
    return(data.table(
      exp   = exp_trait, out = NA_character_,
      exp_i = i, out_i = NA_integer_,
      beta  = NA_real_, se = NA_real_, p = NA_real_,
      nsnp  = n_iv, fdr = NA_real_,
      Notes = sprintf("[Info] Exposure '%s': #IVs ≤ 10 (LOD ≥ 4); skipped.", exp_trait)
    ))
  }
  
  # ---- Normal MR loop for valid exposures ----
  rows <- list(); k <- 1L
  for (j in seq_along(traits)) {
    if (j == i) next
    out_trait <- traits[j]
    
    res <- tryCatch(
      run_mr_from_fastloco(
        exposure_trait = exp_trait,
        outcome_trait  = out_trait,
        res_fast_all   = res_fast_all,
        markers        = markers_iv,
        figure         = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(res)) next
    
    rows[[k]] <- data.table(
      exp   = exp_trait,
      out   = out_trait,
      exp_i = i,
      out_i = j,
      beta  = res$beta_mr,
      se    = res$se_mr,
      p     = res$p_mr,
      nsnp  = res$n_snp,
      Notes = NA_character_
    )
    k <- k + 1L
  }
  
  if (length(rows)) {
    rbindlist(rows, fill = TRUE)
  } else {
    data.table(
      exp   = exp_trait, out = NA_character_,
      exp_i = i, out_i = NA_integer_,
      beta  = NA_real_, se = NA_real_, p = NA_real_,
      nsnp  = n_iv, fdr = NA_real_,
      Notes = sprintf("[Info] Exposure '%s': valid IVs but no MR pairs returned.", exp_trait)
    )
  }
}, mc.cores = cores)

mr_pairs <- rbindlist(res_list, fill = TRUE)

# Global FDR for non-NA p values
if (nrow(mr_pairs) && any(!is.na(mr_pairs$p))) {
  mr_pairs[, fdr := p.adjust(p, method = "fdr")]
}

setorder(mr_pairs, fdr, p)
head(mr_pairs[, .(exp, out, beta, se, p, fdr, nsnp, Notes)])
write.csv(mr_pairs, file = "res_mr_pairs.csv", row.names = F)

mr_pairs_norm <- normalize_mr(mr_pairs, phen)

head(
  mr_pairs_norm %>%
    dplyr::select(exp, out, beta, se, beta_norm, se_norm, p, fdr)
)
write.csv(mr_pairs_norm, file = "res_mr_pairs_norm.csv", row.names = F)