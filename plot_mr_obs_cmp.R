source("./FUN.R")
# load("./res_lmm.rda")
# datadir <- "/Users/xd14188/Desktop/UoB/miceMR/working/data/Farber_Lab_DO_genotype"
# # map <- fread(file.path(datadir, "Marker-positions/GM_markers_pmap.csv"))
# load(file.path(datadir, "QTL2_RData/DO_800_QTL_scan.RData"))
# # load(file.path(datadir, "QTL2_RData/DO_800_allele_probs.RData"))
# load(file.path(datadir, "QTL2_RData/DO_800_kinship_loco.RData"))
# load(file.path(datadir, "QTL2_RData/cross_basic_cleaned.Rdata"))
# 
# # Get the phenotype and the covariate data from the cross object
# # Covariate data
# DO_800_covar_data <- as.data.frame(cross_basic$covar)
# # Add the "Mouse.ID" column in the covar data
# DO_800_covar_data <- DO_800_covar_data %>% rownames_to_column(var="Mouse.ID")
# # Remove extra columns from the covar data
# DO_800_covar_data <- DO_800_covar_data[, c(1:7)]
# # Phenotype data
# DO_800_pheno_data <- as.data.frame(cross_basic$pheno)
# # Encode sex as numeric values in the covariate data
# DO_800_covar_data[,"sex"] = (DO_800_covar_data[,"sex"] == "M")
# # Create a numeric matrix with the covariate data
# DO_800_covar_data_matrix = apply(DO_800_covar_data,2,as.numeric)
# # Make sure the rownames of the numeric matrix corresponds to mouse Ids
# rownames(DO_800_covar_data_matrix) <- rownames(cross_basic$covar)
# # Remove the Mouse.ID and the sac_date columns from the covariate
# DO_800_covar_data_matrix <- DO_800_covar_data_matrix[,- c(1:2)]
# datadir <- "/Users/xd14188/Desktop/UoB/miceMR/working/data/Farber_Lab_DO_genotype"
# geno <- fread(file.path(datadir, "Genotype/DO_800_geno_full.csv"), header = T)
# phen <- fread(file.path(datadir, "Phenotype/DO_800_pheno_full_with_covars.csv"))
# 
# # ---- Preprocessing Genotype ----
# markers <- rownames(DO_800_QTL_scan[[1]])
# geno <- geno[match(markers, geno$Marker), ]
# stopifnot(all(geno$Marker == markers))
# 
# g <- as.matrix(geno[, -1])
# rownames(g) <- geno$Marker
# 
# # Align phenotype and genotype samples
# common_ids <- intersect(phen$Mouse.ID, colnames(g))
# phen <- phen[phen$Mouse.ID %in% common_ids, ]
# g <- g[, match(phen$Mouse.ID, colnames(g))]
# stopifnot(all(phen$Mouse.ID == colnames(g)))
# 
# miceOBS <-  read.csv("/Users/xd14188/Desktop/UoB/miceMR/working/scripts/mr/regOutput/pairwise_reg_results.csv", header = T)
# miceOBS_norm <- normalize_mr_cols(mr_tbl = miceOBS, phen,beta_col = "beta", se_col = "se", p_col = "p", "_norm", overwrite_fdr = F)
# head(miceOBS_norm); dim(miceOBS_norm)
# write.csv(miceOBS_norm, "./res_OBS_norm.csv", row.names = F)

# --- Read data (as before) ---
miceMR  <- read.csv("./res_mr_SF_norm.csv", header = TRUE, check.names = FALSE)
miceOBS <- read.csv("./res_OBS_norm.csv",  header = TRUE, check.names = FALSE)

# --- Keep needed columns and rename for merge ---
miceMR_sub  <- miceMR[, c("exp", "out", "beta_post_norm")]
colnames(miceMR_sub)[3] <- "beta_mr"

miceOBS_sub <- miceOBS[, c("exp", "out", "beta_norm",
                           intersect(c("fdr_post_recalc","fdr_recalc","fdr"), names(miceOBS))[1])]
colnames(miceOBS_sub)[3:4] <- c("beta_obs", "fdr_use")

# Guard against zeros/NA in FDR and compute -log10(FDR)
# Use machine minimum to avoid Inf, then optionally cap extremely large values for plotting
miceOBS_sub$neglog10_fdr <- -log10(pmax(miceOBS_sub$fdr_use, .Machine$double.xmin))
# Soft-cap at the 98th percentile for aesthetics (avoids a few points dominating the size scale)
cap <- stats::quantile(miceOBS_sub$neglog10_fdr, 0.98, na.rm = TRUE)
miceOBS_sub$neglog10_fdr_cap <- pmin(miceOBS_sub$neglog10_fdr, cap)

# --- Merge and clean ---
merged <- merge(miceMR_sub, miceOBS_sub[, c("exp","out","beta_obs","neglog10_fdr_cap")],
                by = c("exp","out"))
merged <- subset(merged, is.finite(beta_mr) & is.finite(beta_obs))

# add some filteration that ctrl the most fdr pairs
merged_ord <- merged[order(merged$neglog10_fdr_cap, decreasing = T), ]
# head(merged_ord, 100)
merged_flt <- merged[which(merged$neglog10_fdr_cap < 100), ]
# --- Plot (Nature Genetics-ready look) ---
# --- Fit lm without intercept: MR effect ~ OBS effect ---
fit <- lm(beta_mr ~ 0 + beta_obs, data = merged_flt)
slope_lm <- coef(fit)[["beta_obs"]]
# slope_lm
# Pick the 10 strongest (largest -log10 FDR) points
top10 <- merged_flt %>%
  arrange(desc(neglog10_fdr_cap)) %>%
  slice(31:50)
top10$label <- paste0(top10$exp, " \u2192 ", top10$out)  # Unicode arrow →

#  Create ggplot
p <- ggplot(merged_flt, aes(x = beta_obs, y = beta_mr)) +
  geom_point(aes(size = neglog10_fdr_cap), alpha = 0.6, shape = 16, stroke = 0) +
  
  # identity line
  geom_abline(slope = 1, intercept = 0, linewidth = 0.6, color = "red") +
  
  # lm line
  geom_abline(slope = slope_lm, intercept = 0,
              linewidth = 0.6, color = "blue", linetype = "dashed") +
  
  # --- ADD LABELS FOR TOP 10 ---
  geom_label_repel(
    data = top10,
    aes(label = label),
    size = 2.2,              # small label text
    max.overlaps = Inf,
    label.size = 0.1,
    fill = "white",
    segment.color = "gray40",
    box.padding = 0.25,
    point.padding = 0.2
  ) +
  
  labs(
    x = expression(paste("Observational effect (", beta[OBS], " = beta_norm)")),
    y = expression(paste("MR effect (", beta[IV], " = beta_post_norm)")),
    size = expression(-log[10]("FDR")),
    title = "Mendelian Randomization vs Observational Regression"
  ) +
  scale_size_continuous(range = c(1.2, 6),
                        breaks = c(1, 2, 3, 5, 10),
                        limits = c(0, cap)) +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width  = unit(0.6, "cm")
  )

print(p)
# --- Optional: save high-res for manuscript ---
# ggsave("MR_vs_OBS_scatter.png", p, width = 85, height = 85, units = "mm", dpi = 600)
# ggsave("MR_vs_OBS_scatter.pdf", p, width = 85, height = 85, units = "mm", useDingbats = FALSE)