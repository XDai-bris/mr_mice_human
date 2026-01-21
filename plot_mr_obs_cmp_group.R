###############################################
## Figure: MR vs Observational regression
## Publication-ready (Nature Genetics style)
##
## Pipeline overview
## -----------------
## Input:
##   - res_mr_SF_norm.csv  : MR results (per exp–out pair)
##   - res_OBS_norm.csv    : Observational results (per exp–out pair)
##
## Steps:
##   1) Create biologically meaningful labels (exp_label, out_label)
##      using a custom grouping function `categorize()`.
##   2) For each label pair (exp_label, out_label), use **MR FDR**
##      to choose *one* representative trait pair (exp, out).
##      --> MR is the "decision maker" for which trait represents the group.
##   3) For that same (exp, out) pair, attach the corresponding
##      observational effect and FDR.
##   4) Prepare a scatter plot:
##         x = observational effect (beta_obs)
##         y = MR effect (beta_mr)
##         point size = -log10(FDR_obs)
##         point color = sign concordance (concordant vs discordant)
##      with identity line, MR~OBS regression line, and
##      labels for key exposure–outcome label pairs.
###############################################

## load functions
source("./FUN.R")

###############################################
## 0) Read full MR & OBS data
###############################################

# Paths to input files (adjust if needed)
mr_file  <- "res_mr_SF_norm.csv"
obs_file <- "res_OBS_norm.csv"

miceMR  <- read.csv(mr_file,  header = TRUE, check.names = FALSE)
miceOBS <- read.csv(obs_file, header = TRUE, check.names = FALSE)

###############################################
## 1) Grouping function: collapse raw traits
##    into biologically meaningful label groups
##
## - Keeps some traits as unique:
##     * body_weight, body_length
## - Groups multiple fat pads together:
##     * RFP, GFP, BFP, FFP -> "fat_pad"
## - Groups femoral geometry traits:
##     * FL, ML, AP -> "femoral_geometry"
## - Groups muscle weights:
##     * soleus_weight, gastroc_weight -> "muscle_weight"
## - For all other traits with a prefix (xxx_yyy):
##     * uses the prefix before "_" as the group name:
##          uCT_Tb.Sp        -> "uCT"
##          bending_max_load -> "bending"
##          histo_TbN        -> "histo"
##
## These labels are used to talk about broad
## exposure/outcome categories in the plot,
## while the actual MR/OBS are still based on
## the underlying specific trait (exp, out).
###############################################

categorize <- function(x) {
  # Explicit groups
  fatpads        <- c("RFP", "GFP", "BFP", "FFP")
  femoral_geom   <- c("FL", "ML", "AP")
  muscle_weights <- c("soleus_weight", "gastroc_weight")
  
  # Traits that should remain their own group,
  # not merged into anything else
  special_unique <- c("body_weight", "body_length")
  
  ifelse(
    x %in% special_unique,
    x,  # its own label
    ifelse(
      x %in% fatpads,        "fat_pad",
      ifelse(
        x %in% femoral_geom, "femoral_geometry",
        ifelse(
          x %in% muscle_weights, "muscle_weight",
          ifelse(
            # default rule: for names like "uCT_Tb.Sp"
            # or "histo_TbN", use the prefix before "_"
            grepl("_", x),
            sub("_.*", "", x),
            # if no "_" and not in any special list,
            # use the original name as label
            x
          )
        )
      )
    )
  )
}

###############################################
## 2) MR side: use MR FDR to choose one trait
##    pair (exp, out) per label group
##
## MR table: miceMR
##  - Columns include: exp, out, beta_post_norm, fdr_post*_*
##
## For each row:
##   exp_label = categorize(exp)
##   out_label = categorize(out)
##
## For each (exp_label, out_label) group, there can be
## multiple specific trait pairs (different exp/out).
##
## We use the **MR FDR** to select ONE representative
## (exp, out) pair per group:
##
##   * OPTION A (default): pick the most significant
##     MR association (smallest MR FDR).
##   * OPTION B (commented below): pick the trait pair
##     whose MR FDR is closest to the median FDR in that group.
##
## This makes MR the "decision maker": MR chooses
## which exact trait pair represents each exposure–
## outcome label combination.
###############################################

## Choose MR FDR column in priority order
if ("fdr_post_recalc" %in% names(miceMR)) {
  fdr_mr_col <- "fdr_post_recalc"
} else if ("fdr_post" %in% names(miceMR)) {
  fdr_mr_col <- "fdr_post"
} else if ("fdr" %in% names(miceMR)) {
  fdr_mr_col <- "fdr"
} else {
  stop("No MR FDR column found in miceMR")
}

# Add labels and MR FDR to MR table
miceMR <- miceMR %>%
  mutate(
    exp_label = categorize(exp),
    out_label = categorize(out),
    fdr_mr    = .data[[fdr_mr_col]]
  )

# Treat missing FDR as non-significant
miceMR$fdr_mr[is.na(miceMR$fdr_mr)] <- 1

# ---- OPTION A: MR chooses the most significant trait per group ----
mr_rep <- miceMR %>%
  group_by(exp_label, out_label) %>%
  # one row per (exp_label, out_label): smallest MR FDR
  slice_min(fdr_mr, with_ties = FALSE) %>%
  ungroup()

# ---- OPTION B (alternative): MR chooses the median trait per group ----
# mr_rep <- miceMR %>%
#   group_by(exp_label, out_label) %>%
#   mutate(fdr_mr_med = median(fdr_mr, na.rm = TRUE)) %>%
#   # pick the row whose MR FDR is closest to the median
#   slice_min(abs(fdr_mr - fdr_mr_med), with_ties = FALSE) %>%
#   ungroup()

# Keep only the columns we need going forward:
# - exp, out: the *specific* trait pair chosen by MR
# - exp_label, out_label: the broad group labels
# - beta_mr, fdr_mr: MR effect size and FDR for that pair
mr_rep_sub <- mr_rep %>%
  transmute(
    exp,
    out,
    exp_label,
    out_label,
    beta_mr = beta_post_norm,
    se_mr   = se_post_norm,
    fdr_mr  = fdr_mr
  )

###############################################
## 3) OBS side: attach OBS FDR & beta to the
##    same MR-chosen trait pairs
##
## We do NOT decide groups based on OBS.
## Instead, we respect MR's choice of the specific
## (exp, out) trait pair within each label group.
##
## Steps:
##   1) From res_OBS_norm.csv, keep beta_norm and FDR
##      for each (exp, out) pair.
##   2) Merge with the MR representative table mr_rep_sub
##      by (exp, out).
##   3) After merging, each row corresponds to:
##        * one exposure–outcome label pair
##        * ONE specific trait pair chosen by MR
##        * MR beta & FDR (beta_mr, fdr_mr)
##        * OBS beta & FDR (beta_obs, fdr_obs)
##   4) For plotting, we use:
##        x = beta_obs
##        y = beta_mr
##        size = -log10(FDR_obs)
##      and labels from (exp_label → out_label).
###############################################

## Choose OBS FDR column in priority order
if ("fdr_recalc" %in% names(miceOBS)) {
  fdr_obs_col <- "fdr_recalc"
} else if ("fdr" %in% names(miceOBS)) {
  fdr_obs_col <- "fdr"
} else {
  stop("No OBS FDR column found in miceOBS")
}

# Keep only the observational effect and FDR
miceOBS_sub <- miceOBS %>%
  transmute(
    exp,
    out,
    beta_obs = beta_norm,
    se_obs   = se_norm,
    fdr_obs  = .data[[fdr_obs_col]]
  )

# Merge: for each MR-chosen (exp, out),
# attach the corresponding OBS result
merged <- mr_rep_sub %>%
  left_join(miceOBS_sub, by = c("exp", "out"))

# Handle missing OBS FDR (treat as non-significant)
merged$fdr_obs[is.na(merged$fdr_obs)] <- 1

# Remove rows with non-finite effect estimates
merged <- merged %>%
  filter(is.finite(beta_mr), is.finite(beta_obs))

# Compute -log10(FDR_mr) and cap extreme values for plotting aesthetics
merged <- merged %>%
  mutate(
    neglog10_fdr     = -log10(pmax(fdr_mr, .Machine$double.xmin)),
    # cap at 98th percentile to avoid a few points dominating
    neglog10_fdr_cap = {
      cap <- stats::quantile(neglog10_fdr, 0.98, na.rm = TRUE)
      pmin(neglog10_fdr, cap)
    }
  )

# This final 'merged' table is what we use for:
#   - regression of MR vs OBS
#   - sign-concordance classification
#   - the publication plot.
merged_flt <- merged

###############################################
# Optional: drop within-label pairs such as
# MAT -> MAT or histo -> histo (if desired)
#
# Uncomment the line below if you only want
# cross-category effects in the figure:

merged_flt <- merged_flt %>%
  filter(exp_label != out_label)
###############################################

###############################################
## 4) Regression and summary statistics
###############################################

# Fit regression: MR effect ~ OBS effect (no intercept)
fit <- lm(beta_mr ~ 0 + beta_obs, data = merged_flt)
slope_lm <- coef(fit)[["beta_obs"]]

# Pearson correlation between MR and OBS effects
r_val <- cor(merged_flt$beta_obs, merged_flt$beta_mr, use = "complete.obs")

###############################################
## 5) Q-test: is beta_MR significantly
##    different from beta_OBS?
##
## Q = (beta_mr - beta_obs)^2 / (se_mr^2 + se_obs^2)
## df = 1  --> p-value from chi-square(1)
###############################################

merged_flt <- merged_flt %>%
  mutate(
    Q     = (beta_mr - beta_obs)^2 / (se_mr^2 + se_obs^2),
    Q_p   = pchisq(Q, df = 1, lower.tail = FALSE),
    Q_grp = dplyr::case_when(
      is.na(Q_p)           ~ "Uncertain SE",
      Q_p < 0.05           ~ "Significantly different (Q < 0.05)",
      TRUE                 ~ "Not significantly different"
    )
  )

merged_flt$Q_grp <- factor(
  merged_flt$Q_grp,
  levels = c("Significantly different (Q < 0.05)",
             "Not significantly different",
             "Uncertain SE")
)
###############################################
## 6) Choose points to label in the figure
##
## We label the top N groups ranked by
## -log10(FDR_obs), showing them as:
##   exp_label → out_label
###############################################

topN <- 10  # how many points to label

label_df <- merged_flt %>%
  arrange(desc(neglog10_fdr_cap)) %>%
  slice(1:topN) %>%
  mutate(
    label = paste0(exp_label, " \u2192 ", out_label)
  )

# Annotation position (top-left-ish in data space)
x_annot <- min(merged_flt$beta_obs, na.rm = TRUE)
y_annot <- max(merged_flt$beta_mr, na.rm = TRUE)

###############################################
## 7) Publication-ready plot:
##    MR vs OBS with sign concordance
##
## - x-axis: observational effect (beta_obs)
## - y-axis: MR effect (beta_mr)
## - point color: sign concordance
## - point size: -log10(FDR_obs)
## - identity line: MR = OBS
## - dotted line: regression MR ~ OBS
## - dashed grey lines at x = 0 and y = 0
## - labels for top N exposure–outcome label pairs
###############################################

p <- ggplot(merged_flt, aes(x = beta_obs, y = beta_mr)) +
  # vertical error bars for MR
  geom_errorbar(
    aes(ymin = beta_mr - 1.96 * se_mr,
        ymax = beta_mr + 1.96 * se_mr,
        color = Q_grp),
    width = 0,
    alpha = 0.4,
    linewidth = 0.3
  ) +
  # horizontal error bars for OBS
  geom_errorbar(
    aes(xmin = beta_obs - 1.96 * se_obs,
        xmax = beta_obs + 1.96 * se_obs,
        color = Q_grp),
    orientation = "y",  # tells ggplot these are horizontal
    width = 0,
    alpha = 0.4,
    linewidth = 0.3
  ) +
  # --- existing points ---
  geom_point(
    aes(size = neglog10_fdr_cap, color = Q_grp),
    alpha = 0.7,
    shape = 16,
    stroke = 0
  ) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.4) +
  # build a tiny data frame so we can map linetype and get a legend
  geom_abline(
    data = data.frame(
      slope     = c(1, slope_lm),
      intercept = c(0, 0),
      line_type = c("y = x", "MR ~ OBS regression")
    ),
    aes(slope = slope, intercept = intercept, linetype = line_type),
    linewidth = 0.6,
    color = "black",
    show.legend = TRUE
  ) +
  
  ## 👇 CHANGED: use label_df instead of top_lab
  geom_label_repel(
    data = label_df,
    aes(x = beta_obs, y = beta_mr, label = label),
    inherit.aes = FALSE,
    size = 2.2,
    max.overlaps = Inf,
    label.size = 0.1,
    fill = "white",
    segment.color = "grey40",
    box.padding = 0.25,
    point.padding = 0.2,
    min.segment.length = 0
  ) +
  
  annotate(
    "text",
    x = x_annot,
    y = y_annot,
    hjust = 0,
    vjust = 1,
    label = paste0("Slope = ", sprintf("%.2f", slope_lm),
                   "\nR = ", sprintf("%.2f", r_val)),
    size = 2.4
  ) +
  labs(
    x = expression(paste("Observational effect (", beta[OBS], ")")),
    y = expression(paste("MR effect (", beta[IV], ")")),
    size = expression(-log[10]("FDR")[MR]),
    color = "Q-test",
    title = NULL
  ) +
  
  # Size scale for FDR
  scale_size_continuous(
    range  = c(1.2, 5.0),
    breaks = c(2, 4, 6, 8, 10)
  ) +
  
  # Simple, interpretable colour scheme:
  #  - greys for concordant
  #  - red for discordant
  scale_color_manual(
    values = c(
      "Significantly different (Q < 0.05)" = "#D55E00",  # reddish
      "Not significantly different"        = "grey30",
      "Uncertain SE"                       = "grey80"
    ),
    drop = TRUE
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("y = x" = "solid",
               "MR ~ OBS regression" = "dotted")
  ) +
  
  coord_equal(expand = TRUE) +
  
  theme_classic(base_size = 7) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title   = element_text(face = "plain"),
    axis.text    = element_text(color = "black"),
    legend.position = "right",
    legend.box      = "vertical",
    legend.title    = element_text(face = "plain"),
    legend.text     = element_text(size = 6),
    axis.ticks      = element_line(linewidth = 0.3),
    plot.margin     = margin(2, 2, 2, 2, "mm")
  ) +
  guides(
    size  = guide_legend(order = 1, override.aes = list(alpha = 1)),
    color = guide_legend(order = 2, override.aes = list(size = 3))
  )

print(p)

###############################################
## 8) Export – single-column Nature Genetics size
###############################################

ggsave(
  "Fig_MR_vs_OBS_concordance_short.pdf",
  p,
  width = 200, height = 150,
  units = "mm",
  useDingbats = FALSE
)

# ggsave(
#   "Fig_MR_vs_OBS_concordance.png",
#   p,
#   width = 85, height = 85,
#   units = "mm",
#   dpi = 300
# )