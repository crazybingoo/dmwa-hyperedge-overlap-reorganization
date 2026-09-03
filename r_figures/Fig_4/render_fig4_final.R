#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
  library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
project_root <- Sys.getenv("PAPER1_PROJECT_ROOT", unset = repo_root)
# The derived motif tables retain their legacy figures_data/Fig_5 location.
fig_source_dir <- file.path(project_root, "figures_data", "Fig_5")
out_dir <- file.path(project_root, "nature_fig", "Fig_4")
dir.create(fig_source_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(repo_root, "r_figures", "shared", "paper1_figure_style_config.R"))

stage_levels <- c("pre-ictal", "early", "mid", "late", "post-ictal")
stage_labels_short <- c(
  "pre-ictal" = "Pre",
  "early" = "Early",
  "mid" = "Mid",
  "late" = "Late",
  "post-ictal" = "Post"
)

target_stage_levels <- c("early", "mid", "late", "post-ictal")
target_stage_labels <- unname(stage_labels_short[target_stage_levels])

motif_vars <- c(
  simple_fraction_mean = "simple motifs",
  reciprocal_fraction_mean = "broad reciprocal family",
  dense_fraction_mean = "dense reciprocal subset"
)

motif_levels <- unname(motif_vars)
motif_palette <- c(
  "simple motifs" = sig("blue"),
  "broad reciprocal family" = sig("red"),
  "dense reciprocal subset" = sig("purple")
)
motif_short_labels <- c(
  "simple motifs" = "simple",
  "broad reciprocal family" = "broad reciprocal",
  "dense reciprocal subset" = "dense reciprocal subset"
)

theme_set(theme_paper1())
panel_tag_theme <- theme(
  plot.tag = element_text(
    size = FIG_PANEL_PT,
    face = "bold",
    colour = paper1_palette$text_dark
  )
)

stage_factor <- function(x) {
  factor(x, levels = stage_levels)
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

bootstrap_mean_ci <- function(x, n_boot = 5000) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(c(ci_low = NA_real_, ci_high = NA_real_))
  boot <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  stats::quantile(boot, c(0.025, 0.975), names = FALSE, na.rm = TRUE) |>
    setNames(c("ci_low", "ci_high"))
}

save_pub_r <- function(plot, fig_id = "Fig_4", width_mm = 183, height_mm = 140, dpi = 600) {
  base <- file.path(out_dir, fig_id)
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white")
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(base, ".tiff"),
    width = w,
    height = h,
    units = "in",
    res = dpi,
    background = "white",
    compression = "lzw"
  )
  print(plot)
  dev.off()
}

stage_path <- file.path(project_root, "stats", "phase13_elbow_motif_stage_summary.csv")
fig5_stage <- read.csv(stage_path, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(phase5 = stage_factor(phase5))

write.csv(
  fig5_stage,
  file.path(fig_source_dir, "fig5_motif_stage_summary_elbow.csv"),
  row.names = FALSE
)

patient_stage <- fig5_stage %>%
  group_by(patient_id, phase5) %>%
  summarise(
    n_cases = n(),
    simple_fraction_mean = mean(simple_fraction_mean, na.rm = TRUE),
    reciprocal_fraction_mean = mean(reciprocal_fraction_mean, na.rm = TRUE),
    dense_fraction_mean = mean(dense_fraction_mean, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  patient_stage,
  file.path(fig_source_dir, "fig5_patient_stage_summary_elbow.csv"),
  row.names = FALSE
)

patient_long <- patient_stage %>%
  pivot_longer(
    cols = all_of(names(motif_vars)),
    names_to = "metric",
    values_to = "fraction"
  ) %>%
  mutate(
    motif = factor(motif_vars[metric], levels = motif_levels),
    phase_label = factor(stage_labels_short[as.character(phase5)], levels = unname(stage_labels_short))
  )

trend_summary <- patient_long %>%
  group_by(phase5, phase_label, motif) %>%
  summarise(
    n_patient = n(),
    mean_fraction = mean(fraction, na.rm = TRUE),
    sem_fraction = sem(fraction),
    .groups = "drop"
  )

write.csv(
  trend_summary,
  file.path(fig_source_dir, "fig5_patient_motif_trend_summary_elbow.csv"),
  row.names = FALSE
)

trend_summary <- trend_summary %>%
  mutate(motif_short = factor(motif_short_labels[as.character(motif)], levels = unname(motif_short_labels)))

pre_ref <- patient_stage %>%
  filter(phase5 == "pre-ictal") %>%
  select(patient_id, all_of(names(motif_vars))) %>%
  rename_with(~ paste0("pre_", .x), all_of(names(motif_vars)))

patient_delta <- patient_stage %>%
  filter(phase5 != "pre-ictal") %>%
  select(patient_id, phase5, all_of(names(motif_vars))) %>%
  inner_join(pre_ref, by = "patient_id") %>%
  mutate(
    simple_delta = simple_fraction_mean - pre_simple_fraction_mean,
    reciprocal_delta = reciprocal_fraction_mean - pre_reciprocal_fraction_mean,
    dense_delta = dense_fraction_mean - pre_dense_fraction_mean
  ) %>%
  select(patient_id, phase5, simple_delta, reciprocal_delta, dense_delta)

write.csv(
  patient_delta,
  file.path(fig_source_dir, "fig5_patient_preictal_deltas_elbow.csv"),
  row.names = FALSE
)

delta_long <- patient_delta %>%
  pivot_longer(
    cols = ends_with("_delta"),
    names_to = "metric",
    values_to = "delta_fraction"
  ) %>%
  mutate(
    metric = sub("_delta$", "_fraction_mean", metric),
    motif = factor(motif_vars[metric], levels = motif_levels),
    phase5 = factor(phase5, levels = target_stage_levels),
    phase_label = factor(stage_labels_short[as.character(phase5)], levels = target_stage_labels)
  )

set.seed(260508)
delta_summary <- delta_long %>%
  group_by(phase5, phase_label, motif) %>%
  group_modify(~ {
    x <- .x$delta_fraction[is.finite(.x$delta_fraction)]
    ci <- bootstrap_mean_ci(x, n_boot = 5000)
    wt <- suppressWarnings(stats::wilcox.test(
      x,
      mu = 0,
      alternative = "two.sided",
      exact = TRUE
    ))
    tibble(
      n_patient = length(x),
      mean_delta = mean(x),
      median_delta = median(x),
      positive_fraction = mean(x > 0),
      ci_low = unname(ci[["ci_low"]]),
      ci_high = unname(ci[["ci_high"]]),
      wilcoxon_W = unname(wt$statistic),
      wilcoxon_exact_p = wt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    bh_adjusted_p = p.adjust(wilcoxon_exact_p, method = "BH"),
    bh_family_n = n(),
    bootstrap_replicates = 5000L,
    bootstrap_seed = 260508L,
    ci_definition = "unadjusted 95% percentile bootstrap CI for the patient-level mean"
  )

write.csv(
  delta_summary,
  file.path(fig_source_dir, "fig5_patient_preictal_delta_summary_elbow.csv"),
  row.names = FALSE
)

delta_axis_limit <- ceiling(max(abs(delta_long$delta_fraction), na.rm = TRUE) * 10) / 10

reciprocal_star_x <- 0.46
reciprocal_star_mean_gap <- reciprocal_star_x - mean(
  delta_summary$mean_delta[as.character(delta_summary$motif) == "broad reciprocal family"]
)
simple_star_x <- mean(
  delta_summary$mean_delta[as.character(delta_summary$motif) == "simple motifs"]
) + reciprocal_star_mean_gap

significance_summary <- delta_summary %>%
  filter(bh_adjusted_p < 0.05) %>%
  mutate(
    significance_label = "*",
    significance_x = ifelse(
      as.character(motif) == "simple motifs",
      simple_star_x,
      reciprocal_star_x
    ),
    significance_hjust = 0.5
  )

# Retain the composition-control calculation as an internal result only. It is
# intentionally not plotted, included in Fig. 4 source data, or promoted to SI.
pre_control <- fig5_stage %>%
  filter(phase5 == "pre-ictal") %>%
  transmute(
    case_id,
    pre_reciprocal_expected = reciprocal_count_mean / all_motif_count_mean,
    pre_dense_expected = dense_count_mean / all_motif_count_mean
  )

control_case <- fig5_stage %>%
  filter(phase5 %in% target_stage_levels) %>%
  left_join(pre_control, by = "case_id") %>%
  transmute(
    case_id,
    target_stage = as.character(phase5),
    reciprocal_empirical = reciprocal_fraction_mean,
    reciprocal_expected = pre_reciprocal_expected,
    dense_empirical = dense_fraction_mean,
    dense_expected = pre_dense_expected
  ) %>%
  pivot_longer(
    cols = c(
      reciprocal_empirical, reciprocal_expected,
      dense_empirical, dense_expected
    ),
    names_to = c("metric_root", ".value"),
    names_pattern = "(reciprocal|dense)_(empirical|expected)"
  ) %>%
  mutate(
    metric = paste0(metric_root, "_fraction"),
    difference = empirical - expected
  ) %>%
  filter(is.finite(empirical), is.finite(expected), is.finite(difference))

control_summary <- control_case %>%
  group_by(target_stage, metric) %>%
  group_modify(~ {
    d <- .x$difference
    wt <- suppressWarnings(stats::wilcox.test(
      d,
      mu = 0,
      alternative = "two.sided",
      exact = TRUE
    ))
    tibble(
      null_model = "elbow_consistent_preictal_composition",
      n_units = length(d),
      empirical_effect = mean(.x$empirical),
      null_mean = mean(.x$expected),
      empirical_minus_null = mean(d),
      median_minus_null = median(d),
      case_percentile_low = unname(quantile(d, 0.025)),
      case_percentile_high = unname(quantile(d, 0.975)),
      positive_fraction = mean(d > 0),
      wilcoxon_W = unname(wt$statistic),
      wilcoxon_exact_p = wt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    bh_adjusted_p = p.adjust(wilcoxon_exact_p, method = "BH"),
    bh_family_n = n(),
    inference_unit = "seizure",
    control_definition = paste(
      "target-stage elbow motif fraction minus the same seizure's",
      "pre-ictal elbow motif-composition expectation"
    )
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      reciprocal_fraction = "broad reciprocal family",
      dense_fraction = "dense reciprocal subset"
    ),
    metric_label = factor(metric_label, levels = c("broad reciprocal family", "dense reciprocal subset")),
    metric_short = factor(
      motif_short_labels[as.character(metric_label)],
      levels = c("broad reciprocal", "dense reciprocal subset")
    ),
    target_stage = factor(target_stage, levels = target_stage_levels),
    stage_label = factor(stage_labels_short[as.character(target_stage)], levels = target_stage_labels),
    stage_num = as.numeric(target_stage),
    metric_offset = ifelse(metric_label == "broad reciprocal family", 0.12, -0.12),
    y_pos = stage_num + metric_offset,
    value_label = sprintf("%+.3f", empirical_minus_null)
  )

control_n <- unique(control_summary$n_units)
if (length(control_n) != 1L) {
  stop("Fig. 4 panel c does not have a consistent eligible-seizure count.")
}

write.csv(
  control_summary,
  file.path(fig_source_dir, "fig5_motif_count_matched_control_elbow.csv"),
  row.names = FALSE
)

write.csv(
  control_case,
  file.path(fig_source_dir, "fig5_motif_count_matched_control_case_elbow.csv"),
  row.names = FALSE
)

motif_trend_plot <- function(motif_name, y_title = FALSE) {
  trend_title <- c(
    "simple motifs" = "simple motifs",
    "broad reciprocal family" = "broad reciprocal family",
    "dense reciprocal subset" = "dense reciprocal subset"
  )[[motif_name]]

  ggplot(
    trend_summary %>% filter(as.character(motif) == motif_name),
    aes(phase_label, mean_fraction, group = 1)
  ) +
    geom_line(colour = motif_palette[[motif_name]], linewidth = 0.54, alpha = 0.92) +
    geom_errorbar(
      aes(ymin = mean_fraction - sem_fraction, ymax = mean_fraction + sem_fraction),
      colour = motif_palette[[motif_name]],
      width = 0.10,
      linewidth = 0.34,
      alpha = 0.88
    ) +
    geom_point(
      colour = motif_palette[[motif_name]],
      shape = 21,
      fill = "white",
      stroke = 0.52,
      size = 1.72
    ) +
    scale_x_discrete(expand = expansion(add = c(0.16, 0.16))) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = c(0.18, 0.18))) +
    labs(
      title = trend_title,
      x = NULL,
      y = if (isTRUE(y_title)) expression(f[F]) else NULL
    ) +
    theme(
      plot.title = element_text(
        size = FIG_TEXT_PT - 0.2,
        face = "bold",
        colour = paper1_palette$text_dark,
        hjust = 0.5,
        margin = margin(0, 0, 2, 0)
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.14, colour = paper1_palette$grid_light),
      axis.title.y = element_text(margin = margin(r = 1.5)),
      plot.margin = margin(1, 5, 5, 4)
    )
}

p5a_simple <- motif_trend_plot("simple motifs", TRUE)
p5a_reciprocal <- motif_trend_plot("broad reciprocal family", FALSE)
p5a_dense <- motif_trend_plot("dense reciprocal subset", FALSE)

p5a_title <- ggplot() +
  labs(title = "(a) Motif-family trajectories") +
  theme_void(base_family = "Arial", base_size = FIG_TEXT_PT) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      size = FIG_PANEL_PT,
      face = "bold",
      colour = paper1_palette$text_dark,
      hjust = 0,
      vjust = 0.5
    ),
    plot.margin = margin(2, 5, 0, 4)
  )

p5b <- ggplot(delta_long, aes(delta_fraction, phase_label, colour = motif)) +
  geom_vline(xintercept = 0, linewidth = 0.34, linetype = "dashed", colour = paper1_palette$annotation_gray) +
  geom_point(
    shape = 21,
    fill = "white",
    stroke = 0.34,
    size = 0.92,
    alpha = 0.34,
    position = position_jitter(width = 0, height = 0.06, seed = 260508)
  ) +
  geom_segment(
    data = delta_summary,
    aes(x = ci_low, xend = ci_high, y = phase_label, yend = phase_label, colour = motif),
    linewidth = 0.46,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = delta_summary,
    aes(mean_delta, phase_label, colour = motif),
    shape = 21,
    fill = "white",
    stroke = 0.60,
    size = 1.85,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = significance_summary,
    aes(
      x = significance_x,
      y = phase_label,
      label = significance_label,
      hjust = significance_hjust
    ),
    colour = paper1_palette$text_dark,
    size = FIG_GEOM_TEXT_SIZE + 0.30,
    fontface = "bold",
    vjust = 0.45,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ motif, nrow = 1) +
  scale_colour_manual(values = motif_palette, guide = "none") +
  scale_x_continuous(
    limits = c(-delta_axis_limit, delta_axis_limit),
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "(b) Paired shifts from pre-ictal baseline",
    x = expression(Delta*f[F]),
    y = NULL
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      size = FIG_PANEL_PT,
      face = "bold",
      colour = paper1_palette$text_dark,
      hjust = 0
    ),
    strip.text = element_text(size = FIG_TEXT_PT, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.16, colour = paper1_palette$grid_light),
    plot.margin = margin(5, 6, 5, 4)
  )

p5c <- ggplot(control_summary, aes(empirical_minus_null, y_pos, colour = metric_label)) +
  geom_vline(xintercept = 0, linewidth = 0.34, linetype = "dashed", colour = paper1_palette$annotation_gray) +
  geom_segment(
    aes(x = 0, xend = empirical_minus_null, y = y_pos, yend = y_pos),
    linewidth = 0.48,
    alpha = 0.88
  ) +
  geom_point(shape = 21, fill = "white", stroke = 0.58, size = 1.95) +
  scale_colour_manual(
    values = motif_palette[c("broad reciprocal family", "dense reciprocal subset")],
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 0.085),
    breaks = seq(0, 0.08, 0.02),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    breaks = seq_along(target_stage_levels),
    labels = target_stage_labels,
    expand = expansion(add = 0.38)
  ) +
  labs(
    title = sprintf("(c) Pre-ictal composition control across %d eligible seizures", control_n),
    x = expression(Delta*M[list(s,F)]^{"emp-null"}),
    y = NULL
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      size = FIG_PANEL_PT,
      face = "bold",
      colour = paper1_palette$text_dark,
      hjust = 0
    ),
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.key.width = unit(9, "pt"),
    legend.spacing.x = unit(4, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, -5, 0),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.14, colour = paper1_palette$grid_light),
    plot.margin = margin(5, 6, 4, 4)
  )

fig5_design <- c(
  area(t = 1, l = 1, b = 1, r = 3),
  area(t = 2, l = 1, b = 2, r = 1),
  area(t = 2, l = 2, b = 2, r = 2),
  area(t = 2, l = 3, b = 2, r = 3),
  area(t = 3, l = 1, b = 3, r = 3)
)

fig5 <- p5a_title + p5a_simple + p5a_reciprocal + p5a_dense + p5b +
  plot_layout(design = fig5_design, heights = c(0.10, 1.00, 1.18)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(fig5, "Fig_4", 183, 126)

source_paths <- c(
  file.path(fig_source_dir, "fig5_motif_stage_summary_elbow.csv"),
  file.path(fig_source_dir, "fig5_patient_stage_summary_elbow.csv"),
  file.path(fig_source_dir, "fig5_patient_motif_trend_summary_elbow.csv"),
  file.path(fig_source_dir, "fig5_patient_preictal_deltas_elbow.csv"),
  file.path(fig_source_dir, "fig5_patient_preictal_delta_summary_elbow.csv")
)

if (!all(file.exists(source_paths))) {
  stop("One or more Fig. 4 source-data files are missing.")
}

combined_source_data <- bind_rows(lapply(source_paths, function(path) {
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  normalized_path <- gsub("\\\\", "/", path)
  normalized_root <- gsub("\\\\", "/", project_root)
  dat$.source_file <- substring(normalized_path, nchar(normalized_root) + 2L)
  dat$.source_row <- seq_len(nrow(dat))
  dat %>% relocate(.source_file, .source_row)
}))

source_data_path <- file.path(out_dir, "Fig_4_source_data.csv")
write.csv(combined_source_data, source_data_path, row.names = FALSE)

figure_paths <- file.path(out_dir, paste0("Fig_4", c(".svg", ".pdf", ".png", ".tiff")))
export_manifest <- data.frame(
  fig_id = "Fig_4",
  file = gsub("\\\\", "/", figure_paths),
  width_mm = 183,
  height_mm = 126,
  dpi = c(NA, NA, 600, 600),
  exists = file.exists(figure_paths),
  size_bytes = file.info(figure_paths)$size,
  stringsAsFactors = FALSE
)
write.csv(export_manifest, file.path(out_dir, "Fig_4_export_manifest.csv"), row.names = FALSE)

if (!all(export_manifest$exists & export_manifest$size_bytes > 0)) {
  stop("One or more Fig. 4 exports are missing or empty.")
}

message("Wrote Fig. 4 redraw to: ", out_dir)
