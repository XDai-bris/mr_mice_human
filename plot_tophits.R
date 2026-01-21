# This src is for plot illustration figure for peaking up tophits with LOD
# output:tophits.png

source("./FUN.R")

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
a  <- as.numeric(DO_800_QTL_scan[[2]])
a1 <- get_tophits(a)       # peak indices
lod_thresh <- 4

# Build a data frame for ggplot
df <- data.frame(
  pos     = seq_along(a),
  lod     = a,
  is_peak = seq_along(a) %in% a1
)

ggplot(df, aes(x = pos, y = lod)) +
  geom_line() +
  geom_hline(yintercept = lod_thresh, linetype = "dashed") +
  geom_point(
    data = subset(df, is_peak),
    color = "red",
    size = 2
  ) +
  labs(
    x = "Position (index)",
    y = "LOD score",
    title = "QTL scan with detected peaks"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )