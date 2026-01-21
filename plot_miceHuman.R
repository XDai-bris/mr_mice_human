summary_tbl <- read.csv("./mice_human_mr_sum_flt.csv", header = T)
## adding the description for the mice phenotype
summary_tbl$mice_exp[which(summary_tbl$mice_exp == "RFP")] <- "RFP(Retro peritoneal fat pad)"
summary_tbl$mice_exp[which(summary_tbl$mice_exp == "FFP")] <- "FFP(Femoral fat pad)"
summary_tbl$mice_exp[which(summary_tbl$mice_exp == "FL")] <- "FL(Full femoral length)"
summary_tbl$mice_exp[which(summary_tbl$mice_exp == "ML")] <- "ML(Medial-lateral femoral width)"

summary_tbl$mice_out[which(summary_tbl$mice_out == "RFP")] <- "RFP(Retro peritoneal fat pad)"
summary_tbl$mice_out[which(summary_tbl$mice_out == "FFP")] <- "FFP(Femoral fat pad)"
summary_tbl$mice_out[which(summary_tbl$mice_out == "FL")] <- "FL(Full femoral length)"
summary_tbl$mice_out[which(summary_tbl$mice_out == "ML")] <- "ML(Medial-lateral femoral width)"

# add tmp filteration to remove FL(Full femoral length) and ML(Medial-lateral femoral width)
summary_tbl_filtered <- summary_tbl %>%
  filter(!(mice_exp == "FL(Full femoral length)" | mice_out == "FL(Full femoral length)"))
summary_tbl <- summary_tbl_filtered
summary_tbl_filtered <- summary_tbl %>%
  filter(!(mice_exp == "ML(Medial-lateral femoral width)" | mice_out == "ML(Medial-lateral femoral width)"))
summary_tbl <- summary_tbl_filtered
## --- 0) Build forest_df with two-line labels (as before) --------------------
forest_df <- bind_rows(
  summary_tbl %>%
    transmute(
      mice_exp, mice_out, human_exp, human_out,
      beta = human_beta, se = human_se, source = "Human", human_p
    ),
  summary_tbl %>%
    transmute(
      mice_exp, mice_out, human_exp, human_out,
      beta = mice_beta, se = mice_se, source = "Mouse", human_p
    )
) %>%
  mutate(
    ci_lo = beta - 1.96 * se,
    ci_hi = beta + 1.96 * se,
    pair_label = paste0(
      mice_exp, " → ", mice_out, "\n",
      "(", human_exp, " → ", human_out, ")"
    )
  ) %>%
  filter(is.finite(beta), is.finite(se))

## --- 1) Order pairs by human_p (smallest first) -----------------------------
pair_order <- summary_tbl %>%
  arrange(human_p) %>%
  mutate(pair_label = paste0(
    mice_exp, " → ", mice_out, "\n",
    "(", human_exp, " → ", human_out, ")"
  )) %>%
  pull(pair_label) %>%
  unique()

forest_df <- forest_df %>%
  mutate(pair_label = factor(pair_label, levels = pair_order))

## --- 2) Split into chunks of 15 ---------------------------------------------
max_per_fig <- 15L
labels_all  <- levels(forest_df$pair_label)
n_pairs     <- length(labels_all)
chunks <- split(labels_all, ceiling(seq_along(labels_all) / max_per_fig))

## --- 3) Plotting helper ------------------------------------------------------
plot_chunk <- function(df_sub, title_suffix = "") {
  ggplot(
    df_sub,
    aes(y = pair_label, x = beta, xmin = ci_lo, xmax = ci_hi, color = source)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_pointrange(position = position_dodge(width = 0.5), linewidth = 0.6) +
    scale_color_manual(values = c(Mouse = "#1f77b4", Human = "#d62728")) + # blue/red
    labs(
      x = "Effect size (β ± 95% CI)",
      y = NULL,
      color = NULL,
      title = paste0("Mouse vs Human MR estimates per pair", title_suffix)
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 8.5, lineheight = 0.9)
    )
}

# --- 4) Generate and save separate PNGs per chunk (fixed) --------------------
for (chunk_idx in seq_along(chunks)) {
  lbls <- chunks[[chunk_idx]]
  df_sub <- forest_df %>%
    dplyr::filter(pair_label %in% lbls) %>%
    droplevels() %>%
    dplyr::mutate(pair_label = factor(pair_label, levels = rev(levels(pair_label))))
  
  start_idx <- (chunk_idx - 1L) * max_per_fig + 1L
  end_idx   <- min(chunk_idx * max_per_fig, n_pairs)
  title_suf <- sprintf(" (pairs %d–%d of %d)", start_idx, end_idx, n_pairs)
  
  p <- plot_chunk(df_sub, title_suffix = title_suf)
  png_file <- sprintf("./mr_forest_pairs_fig/mr_forest_pairs_%03d-%03d.png", start_idx, end_idx)
  ggsave(
    png_file, p,
    width = 8,
    height = 0.35 * nlevels(df_sub$pair_label) + 2,
    dpi = 300
  )
  message("Saved: ", png_file)
}
## --- 5) OPTIONAL: one multi-page PDF instead of multiple PNGs ---------------
# pdf("mr_forest_all_chunks.pdf", width = 8,
#     height = 0.35 * min(max_per_fig, n_pairs) + 2)
# for (idx in seq_along(chunks)) {
#   lbls   <- chunks[[idx]]
#   df_sub <- forest_df %>% filter(pair_label %in% lbls) %>%
#     droplevels() %>%
#     mutate(pair_label = factor(pair_label, levels = rev(levels(pair_label))))
#   start_idx <- (idx - 1L) * max_per_fig + 1L
#   end_idx   <- min(idx * max_per_fig, n_pairs)
#   title_suf <- sprintf(" (pairs %d–%d of %d)", start_idx, end_idx, n_pairs)
#   print(plot_chunk(df_sub, title_suffix = title_suf))
# }
# dev.off()



# ---- Parameters ----
use_numeric_tags <- FALSE 
# ---- Prep ----
df <- summary_tbl %>%
  mutate(
    # 95% CI for error bars
    human_ci = 1.96 * human_se,
    mice_ci  = 1.96 * mice_se,
    # Label from run_name: "exp:body_length_out:body_weight" -> "body_length → body_weight"
    label_full = run_name %>%
      str_remove("^exp:") %>%
      str_replace("_out:", " \u2192 "),
    # Distance from parity line y = x (proportional to perpendicular distance)
    dev_parity = abs(mice_beta - human_beta)
  )

# Pick top 3 farthest points
df_top3 <- df %>%
  arrange(desc(dev_parity)) %>%
  slice_head(n = 15) %>%
  mutate(tag = paste0("[", row_number(), "]"))

# Choose label variable depending on tag mode
if (use_numeric_tags) {
  df <- df %>% left_join(df_top3 %>% select(run_name, tag), by = "run_name")
  df_top3$label_plot <- df_top3$tag
} else {
  df_top3$label_plot <- df_top3$label_full
}

# ---- Stats for annotation ----
fit <- lm(mice_beta ~ human_beta, data = df)
fit_s <- summary(fit)
slope  <- unname(coef(fit)[2])
slope_se <- fit_s$coefficients[2, 2]
r_pearson <- cor(df$human_beta, df$mice_beta, use = "complete.obs")

ann_x <- min(df$human_beta, na.rm = TRUE) + 0.02*diff(range(df$human_beta, na.rm = TRUE))
ann_y <- max(df$mice_beta,  na.rm = TRUE) - 0.02*diff(range(df$mice_beta,  na.rm = TRUE))
ann_text <- sprintf("OLS slope = %.2f \u00B1 %.2f\nPearson r = %.2f", slope, slope_se, r_pearson)

# ---- Plot ----
p <- ggplot(df, aes(x = human_beta, y = mice_beta)) +
  # Parity line y = x
  geom_abline(slope = 1, intercept = 0, linewidth = 0.6, linetype = "22", color = "grey40") +
  
  # Error bars (subtle)
  geom_errorbar(aes(ymin = mice_beta - mice_ci, ymax = mice_beta + mice_ci),
                width = 0, linewidth = 0.45, color = "grey40") +
  geom_segment(aes(x = human_beta - human_ci, xend = human_beta + human_ci,
                   y = mice_beta, yend = mice_beta),
               linewidth = 0.45, color = "grey40") +
  
  # Base points (open circles)
  geom_point(shape = 21, size = 2.8, stroke = 0.8, fill = "white", color = "black") +
  
  # Highlight top-3 (filled)
  geom_point(data = df_top3, shape = 21, size = 3.6, stroke = 1, fill = "black", color = "black") +
  
  # Labels only for top-3, subtle connectors
  ggrepel::geom_text_repel(
    data = df_top3,
    aes(label = label_plot),
    size = 3,
    color = "grey20",
    segment.color = "grey60",
    segment.size = 0.3,
    box.padding = 0.28,
    point.padding = 0.22,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  
  # Linear fit (black)
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "black") +
  
  # Corner annotation (stats)
  annotate("label", x = ann_x, y = ann_y, label = ann_text, hjust = 0, vjust = 1,
           label.size = 0.25, size = 3) +
  
  labs(
    title = "Concordance of human vs. mouse effect sizes",
    x = expression(paste("Human effect size  (", beta, " ± 1.96·SE)")),
    y = expression(paste("Mouse effect size  (", beta, " ± 1.96·SE)"))
  ) +
  coord_fixed(ratio = 1) +
  
  # Minimalist, high-contrast grayscale styling
  theme_classic(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "grey98", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 6)),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    axis.ticks.length = unit(2, "pt")
  )

p
ggsave("human_vs_mouse_betas.pdf", p, width = 180, height = 180, units = "mm", useDingbats = FALSE)