# This scr is for compare different linear model to generate effects
# First using LMM-exact model with simple regression
# Second using LMM-exact model with LMM-fast model
# Input: lmm_effec.R
# Output:First > cmp_ols_lmm.png cmp_ols_lmm_snp.png
#        Second > lmm_cmp.png lmm_cmp_snp.png
# Build long table: two rows per marker (fit1 vs fast)
excl <- NULL # or NULL
ind <- rep(TRUE, length(markers))
effs_lmm <- res_fast
if (length(excl) != 0) { ind[excl] <- FALSE }
print(ind)
df_long <- tibble(
  marker    = markers[ind],
  beta_fit1 = effs_lmm$beta_orig[ind],
  se_fit1   = effs_lmm$se_orig[ind],
  beta_fast = effs$bhat[ind],
  se_fast   = effs$se[ind]
) %>%
  pivot_longer(
    cols = c(beta_fit1, beta_fast, se_fit1, se_fast),
    names_to = c(".value", "method"),
    names_pattern = "(beta|se)_(fit1|fast)"
  ) %>%
  mutate(
    method = recode(method,
                    fit1 = "fit1 (mixed model)",
                    fast = "fast_assoc (OLS)")
  ) %>%
  filter(is.finite(beta), is.finite(se))  # drop incomplete rows

# Position dodge so the two methods sit side-by-side per marker
pd <- position_dodge(width = 0.6)

ggplot(df_long, aes(x = marker, y = beta, color = method)) +
  geom_errorbar(aes(ymin = beta - se, ymax = beta + se),
                width = 0.3, position = pd) +
  geom_point(size = 2.8, position = pd) +
  scale_color_manual(values = c(
    "fit1 (mixed model)" = "#d62728",  # red
    "fast_assoc (OLS)"   = "#1f77b4"   # blue
  )) +
  labs(
    title = "Per-marker effect sizes with SE (fit1 vs fast_assoc)",
    x = "Marker",
    y = expression(hat(beta)),
    color = "Method"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



df_plot <- data.frame(
  marker    = markers[ind],
  beta_fit1 = effs_lmm$beta_orig[ind],
  se_fit1   = effs_lmm$se_orig[ind],
  beta_fast = effs$bhat[ind],
  se_fast   = effs$se[ind]
)

ggplot(df_plot, aes(x = beta_fit1, y = beta_fast)) +
  # Horizontal SE bars for fit1 (red)
  geom_errorbarh(aes(xmin = beta_fit1 - se_fit1, xmax = beta_fit1 + se_fit1),
                 color = "#d62728", alpha = 0.6, height = 0.0) +
  # Vertical SE bars for fast_assoc (blue)
  geom_errorbar(aes(ymin = beta_fast - se_fast, ymax = beta_fast + se_fast),
                color = "#1f77b4", alpha = 0.6, width = 0.0) +
  # Central points
  geom_point(size = 3, color = "black") +
  # Reference line (perfect agreement)
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  # Axis limits
  # coord_cartesian(xlim = c(-40, 40), ylim = c(-40, 40)) +
  labs(
    title = "Comparison of β estimates (fit1 vs fast_assoc)",
    x = expression(beta["fit1"]~"(± SE, red horizontal bars)"),
    y = expression(beta["fast_assoc"]~"(± SE, blue vertical bars)")
  ) +
  theme_bw(base_size = 14)



# ------------------------------------------------------------------------------
# COMPARE LMM-OPTIMAL WITH LMM-FAST
# ------------------------------------------------------------------------------
excl <- NULL # or NULL
ind <- rep(TRUE, length(markers))
if (length(excl) != 0) { ind[excl] <- FALSE }
print(ind)
df_long <- tibble(
  marker    = markers[ind],
  beta_fit1 = effs_lmm$beta_orig[ind],
  se_fit1   = effs_lmm$se_orig[ind],
  beta_fast = res_fast$beta_orig[ind],
  se_fast   = res_fast$se_orig[ind]
) %>%
  pivot_longer(
    cols = c(beta_fit1, beta_fast, se_fit1, se_fast),
    names_to = c(".value", "method"),
    names_pattern = "(beta|se)_(LMM_loco_opt|LMM_fast_gls)"
  ) %>%
  mutate(
    method = recode(method,
                    fit1 = "LMM_loco_opt",
                    fast = "LMM_fast_gls")
  ) %>%
  filter(is.finite(beta), is.finite(se))  # drop incomplete rows

# Position dodge so the two methods sit side-by-side per marker
pd <- position_dodge(width = 0.6)

ggplot(df_long, aes(x = marker, y = beta, color = method)) +
  geom_errorbar(aes(ymin = beta - se, ymax = beta + se),
                width = 0.3, position = pd) +
  geom_point(size = 2.8, position = pd) +
  scale_color_manual(values = c(
    "LMM_loco_opt" = "#d62728",  # red
    "LMM_fast_gls"   = "#1f77b4"   # blue
  )) +
  labs(
    title = "Per-marker effect sizes with SE (LMM_loco_opt vs LMM_fast_gls)",
    x = "Marker",
    y = expression(hat(beta)),
    color = "Method"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

df_plot <- data.frame(
  marker    = markers[ind],
  beta_fit1 = effs_lmm$beta_orig[ind],
  se_fit1   = effs_lmm$se_orig[ind],
  beta_fast = res_fast$beta_orig[ind],
  se_fast   = res_fast$se_orig[ind]
)

ggplot(df_plot, aes(x = beta_fit1, y = beta_fast)) +
  # Horizontal SE bars for fit1 (red)
  geom_errorbarh(aes(xmin = beta_fit1 - se_fit1, xmax = beta_fit1 + se_fit1),
                 color = "#d62728", alpha = 0.6, height = 0.0) +
  # Vertical SE bars for fast_assoc (blue)
  geom_errorbar(aes(ymin = beta_fast - se_fast, ymax = beta_fast + se_fast),
                color = "#1f77b4", alpha = 0.6, width = 0.0) +
  # Central points
  geom_point(size = 3, color = "black") +
  # Reference line (perfect agreement)
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  # Axis limits
  # coord_cartesian(xlim = c(-40, 40), ylim = c(-40, 40)) +
  labs(
    title = "Comparison of β estimates (LMM_loco_opt vs LMM_fast_gls)",
    x = expression(beta["LMM_loco_opt"]~"(± SE, red horizontal bars)"),
    y = expression(beta["LMM_fast_gls"]~"(± SE, blue vertical bars)")
  ) +
  theme_bw(base_size = 14)