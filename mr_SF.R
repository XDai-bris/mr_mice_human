source("./FUN.R")
load("./res_lmm.rda")
# Choose columns: prefer beta_orig/se_orig, else beta/se
.pick_effect_cols <- function(dt) {
  has_orig <- all(c("beta_orig","se_orig") %in% names(dt))
  if (has_orig) {
    dt[, .(id = get(if ("marker" %in% names(dt)) "marker" else "snp"),
           chr = chr, maf = maf,
           bhat = beta_orig, se = se_orig)]
  } else { print("can not find beta_orig and se_orig")
  }
}


system.time({ # random try a pair of MR 
  exp_trait  <- names(phen)[round(runif(1,1,82))]
  markers_iv <- get_markers_for_trait(exp_trait, DO_800_QTL_scan, lod_thresh = 4)
  
  res_sf <- run_mr_with_steiger_markers(
    exp_trait, out_trait = names(phen)[round(runif(1,1,82))],
    res_fast_lmm = res_fast_all,
    markers_iv   = markers_iv,
    N = 800,
    sf_p_threshold = 1,
    .figure = TRUE
  )
  res_sf 
})



#### ALL pairs
# ----- CONFIG -----
traits        <- names(DO_800_QTL_scan)     # all 82 traits
lod_thresh_iv <- 4                          # IV threshold for find_peaks markers
N_sample      <- 800                        # sample size used in Steiger p calc
sf_p_thr      <- 1                          # Steiger p threshold (same as your example)
n_cores       <- 6                          # adjust as you like

# Build all exposure–outcome pairs (exclude self-pairs)
pair_grid <- CJ(exp = traits, out = traits)[exp != out]

# Parallel runner over pairs
pboptions(type = "timer")  # progress timer

mr_sf_list <- pbapply::pblapply(
  X = seq_len(nrow(pair_grid)),
  FUN = function(i) run_one_pair(pair_grid$exp[i], pair_grid$out[i]),
  cl = n_cores  # pbapply will parallelize via cluster
)
# Bind results
mr_sf_all <- rbindlist(mr_sf_list, fill = TRUE)

# Add FDR (post and pre) for convenience
mr_sf_all[, fdr_post := p.adjust(p_post, method = "BH")]
mr_sf_all[, fdr_pre  := p.adjust(p_pre,  method = "BH")]

# Order by outcome significance (post-SF)
setorder(mr_sf_all, fdr_post, p_post, exp, out)

# Peek
print(dim(mr_sf_all))
print(head(mr_sf_all, 10))

# Normalise the beta and se by  sd_exp / sd_out
mr_sf_norm <- normalize_mr_cols(
  mr_tbl       = mr_sf_all,
  phen         = phen,
  beta_col     = "beta_post",
  se_col       = "se_post",
  p_col        = "p_post",
  suffix       = "_norm",
  overwrite_fdr = FALSE
)
write.csv(mr_sf_norm, file = "res_mr_SF_norm.csv", row.names = F)

# ----- Q_pre_post figure -----
# Make pre vs post SF figure
# Prepare data
plot_df <- mr_sf_all %>%
  mutate(
    delta_q   = cochran_Q_post - cochran_Q_pre,
    abs_delta = abs(delta_q),
    size_val  = -log10(cochran_p_pre),  # point size based on pre-filtering significance
    pair_lab  = paste(exp, "→", out)
  )

# Label the biggest movers
n_label <- min(8, nrow(plot_df))
label_df <- plot_df %>%
  arrange(desc(abs_delta)) %>%
  slice_head(n = n_label)

# Correlation for annotation
rho <- cor(plot_df$cochran_Q_pre, plot_df$cochran_Q_post,
           method = "spearman", use = "complete.obs")
rho_lab <- paste0("Spearman \u03C1 = ", sprintf("%.2f", rho))

# Plot
p <- ggplot(plot_df, aes(x = cochran_Q_pre, y = cochran_Q_post)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.6, color = "grey40") +
  geom_point(aes(color = abs_delta, size = size_val), alpha = 0.9) +
  geom_text_repel(
    data = label_df,
    aes(label = pair_lab),
    seed = 42,
    box.padding = 0.5,
    point.padding = 0.25,
    segment.size = 0.25,
    max.overlaps = Inf
  ) +
  coord_equal() +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  scale_color_viridis_c(name = "|ΔQ| (post − pre)", option = "plasma", end = 0.95) +
  scale_size_continuous(name = expression(-log[10](italic(p)[pre])),
                        range = c(3, 10), limits = c(0, NA)) +
  labs(
    title = "Cochran’s Q: Before vs After Steiger Filtering",
    subtitle = "Color = |ΔQ|; Size = -log_10(Q_Pval_pre); dashed line is y = x",
    x = "Cochran’s Q (pre)",
    y = "Cochran’s Q (post)",
    caption = "Labeled points: largest |ΔQ| changes."
  ) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5, label = rho_lab) +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey30"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank()
  ); p
ggsave("cochran_Q_pre_vs_post.png", plot = p, width = 12, height = 6, dpi = 500)



