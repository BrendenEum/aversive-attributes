# Parameter-recovery figure for the Bayesian attribute-level model: recovered
# against true values of the per-subject noise sigma (one point per subject
# and dataset, by study), the prior position q (one point per range condition
# and dataset) and the loss aversion lambda (shared by both range conditions,
# so one point per dataset, "Pooled"). The q and lambda panels are both titled
# "Pooled", centred like the Study 1 / Study 2 strips over the sigma facets.
# Every panel has linear axes with the
# same range and ticks on x and y, and a dashed identity line drawn over the
# points. One legend (Narrow, Wide, Pooled) sits at the right of the figure.
# Used by code/parameter_recovery/recover_bayesian_attributes.R, which saves
# it to paper/figures/parameter_recovery_bayesian.pdf, and shown in
# notebooks/bayesian_model.Rmd.
#   sig_df  columns dataset, Study, Subject, WN, true, recovered
#   grp_df  columns dataset, WN, true_q, rec_q, true_lambda, rec_lambda
# Needs ggplot2 and patchwork attached, and the options in MyPlotOptions.R.
plot_parameter_recovery <- function(sig_df, grp_df) {
  # One key for all three panels: the range condition of a point, or "Pooled"
  # for lambda. The q panel carries the legend, the other two hide theirs, and
  # patchwork collects that one legend to the right of the figure. (With the
  # installed patchwork 1.3.0, `& theme()` on a patchwork does not work, so
  # nothing is set on the patchwork itself.)
  key_levels <- c("Narrow", "Wide", "Pooled")
  key_scales <- list(
    scale_color_manual(name = "Range", values = c(.wn_colors, Pooled = "black"), limits = key_levels, drop = FALSE),
    scale_shape_manual(name = "Range", values = c(.wn_shapes, Pooled = 16),      limits = key_levels, drop = FALSE)
  )
  # The identity line, added after the points so it is drawn on top of them;
  # lighter in the q and lambda panels so their few large markers stand out.
  identity_line <- function(alpha = 1)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = .linewidth, color = "black", alpha = alpha)

  # Noise per subject, by study. Both axes show 0 to 40; the few recovered
  # values above 40 fall outside the view but stay in the data (coord_cartesian
  # zooms, it does not drop points).
  p_sigma <- ggplot(sig_df, aes(x = true, y = recovered, color = WN, shape = WN)) +
    facet_wrap(~ Study, labeller = as_labeller(c("1" = "Study 1", "2" = "Study 2"))) +
    geom_point(alpha = 0.2, size = 1) +
    identity_line() +
    scale_x_continuous(breaks = seq(0, 40, 10)) +
    scale_y_continuous(breaks = seq(0, 40, 10)) +
    key_scales +
    labs(x = expression("True " * sigma), y = expression("Recovered " * sigma)) +
    .myPlot +
    coord_cartesian(xlim = c(0, 40), ylim = c(0, 40)) +
    theme(strip.background = element_blank(), strip.text = element_text(size = .title_pt * .fontscale)) +
    guides(color = "none", shape = "none")

  # Prior position, one point per range condition and dataset.
  lims_q <- range(c(grp_df$true_q, grp_df$rec_q))
  p_q <- ggplot(grp_df, aes(x = true_q, y = rec_q, color = WN, shape = WN)) +
    geom_point(size = 6, show.legend = TRUE) +   # TRUE draws a key for "Pooled" too, which this panel has no points for
    identity_line(alpha = 0.75) +
    scale_x_continuous(limits = lims_q) +
    scale_y_continuous(limits = lims_q) +
    key_scales +
    labs(title = "Pooled", x = "True q", y = "Recovered q") +
    .myPlot +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

  # Loss aversion, shared by both conditions: one point per dataset.
  lam_df <- grp_df[grp_df$WN == "Narrow", ]
  lims_l <- range(c(lam_df$true_lambda, lam_df$rec_lambda))
  p_lambda <- ggplot(lam_df, aes(x = true_lambda, y = rec_lambda, color = "Pooled", shape = "Pooled")) +
    geom_point(size = 6) +
    identity_line(alpha = 0.75) +
    scale_x_continuous(limits = lims_l) +
    scale_y_continuous(limits = lims_l) +
    key_scales +
    labs(title = "Pooled", x = expression("True " * lambda), y = expression("Recovered " * lambda)) +
    .myPlot +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = "none", shape = "none")

  (p_sigma / (p_q | p_lambda)) + plot_layout(guides = "collect")
}
