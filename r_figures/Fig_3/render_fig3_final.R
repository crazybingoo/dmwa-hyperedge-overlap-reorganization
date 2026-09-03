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
stats_root <- file.path(project_root, "stats")
layout_mode <- Sys.getenv("FIG3_LAYOUT", unset = "compact")
if (!layout_mode %in% c("compact", "bwide", "b3col", "b3col_adist", "b3col_adotmap", "b3col_acellmatrix", "b3col_avstrips", "b3col_apairbars")) {
  stop("FIG3_LAYOUT must be 'compact', 'bwide', 'b3col', 'b3col_adist', 'b3col_adotmap', 'b3col_acellmatrix', 'b3col_avstrips' or 'b3col_apairbars'", call. = FALSE)
}
candidate_variant <- if (layout_mode == "bwide") {
  "individual_data_refined_bwide"
} else if (layout_mode == "b3col") {
  "individual_data_refined_b3col"
} else if (layout_mode == "b3col_adist") {
  "individual_data_refined_b3col_adist"
} else if (layout_mode == "b3col_adotmap") {
  "individual_data_refined_b3col_adotmap"
} else if (layout_mode == "b3col_acellmatrix") {
  "individual_data_refined_b3col_acellmatrix"
} else if (layout_mode == "b3col_avstrips") {
  "individual_data_refined_b3col_avstrips"
} else if (layout_mode == "b3col_apairbars") {
  "individual_data_refined_b3col_apairbars"
} else {
  "individual_data_refined"
}
candidate_root <- file.path(
  project_root, "nature_fig", "Fig_3", "candidates", candidate_variant
)
dir.create(candidate_root, showWarnings = FALSE, recursive = TRUE)

source(file.path(repo_root, "r_figures", "shared", "paper1_figure_style_config.R"))

theme_set(
  theme_paper1(base_size = FIG_TEXT_PT) +
    theme(
      plot.title = element_text(
        size = FIG_TEXT_PT,
        face = "bold",
        hjust = 0,
        colour = paper1_palette$text_dark
      )
    )
)

stage_levels <- c("pre-ictal", "early", "mid", "late", "post-ictal")
target_stage_levels <- c("early", "mid", "late", "post-ictal")
stage_labels_short <- c(
  "pre-ictal" = "Pre",
  "early" = "Early",
  "mid" = "Mid",
  "late" = "Late",
  "post-ictal" = "Post"
)
stage_palette_fig3 <- c(
  "pre-ictal" = "#6F8FB7",
  "early" = "#76A973",
  "mid" = "#D49B43",
  "late" = "#C96F6A",
  "post-ictal" = "#9A88BE"
)
stage_outline_palette <- stage_palette_fig3[target_stage_levels]

metric_keys <- c(
  "S",
  "source_sink_gini_gap",
  "sink_minus_source_entropy_gap",
  "sink_minus_source_pr_gap"
)
metric_plain <- c(
  "S" = "S",
  "source_sink_gini_gap" = "D[G]",
  "sink_minus_source_entropy_gap" = "D[H]",
  "sink_minus_source_pr_gap" = "D[PR]"
)
metric_bar_labels <- c(
  "S" = "S",
  "source_sink_gini_gap" = "D[G]",
  "sink_minus_source_entropy_gap" = "D[H]",
  "sink_minus_source_pr_gap" = "D[PR]"
)
line_metric_keys <- metric_keys[-1]
line_metric_colours_short <- c(
  "D[G]" = paper1_palette$primary_teal,
  "D[H]" = "#6F7782",
  "D[PR]" = "#A9823A"
)

stage_factor <- function(x) {
  x <- ifelse(x == "post_ictal", "post-ictal", x)
  factor(x, levels = stage_levels)
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

save_pub_r <- function(plot, base, width_mm = 183, height_mm = 148, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  draw_plot <- function() {
    if (is.function(plot)) plot() else print(plot)
  }

  svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in, bg = "white")
  draw_plot()
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(base, ".pdf"), width = width_in, height = height_in,
    family = "Arial", bg = "white"
  )
  draw_plot()
  grDevices::dev.off()

  ragg::agg_png(
    paste0(base, ".png"), width = width_in, height = height_in,
    units = "in", res = dpi, background = "white"
  )
  draw_plot()
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(base, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, background = "white", compression = "lzw"
  )
  draw_plot()
  grDevices::dev.off()

  data.frame(
    fig_id = "Fig_3_individual_data_refined",
    file = paste0(base, c(".svg", ".pdf", ".png", ".tiff")),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = c(NA, NA, dpi, dpi),
    stringsAsFactors = FALSE
  )
}

patient_stage <- read.csv(
  file.path(stats_root, "phase9_elbow_threshold_patient_stage_features.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(stage = stage_factor(stage))

long_metrics <- patient_stage %>%
  select(patient_id, stage, all_of(metric_keys)) %>%
  pivot_longer(cols = all_of(metric_keys), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = metric_keys),
    metric_label = factor(metric_plain[as.character(metric)], levels = metric_plain[metric_keys])
  )

stage_desc <- long_metrics %>%
  group_by(stage, metric, metric_label) %>%
  summarise(
    n_patients = n_distinct(patient_id),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    sem = sem(value),
    median = median(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

diff_long <- long_metrics %>%
  group_by(patient_id, metric) %>%
  mutate(pre_value = value[stage == "pre-ictal"][1]) %>%
  ungroup() %>%
  filter(stage != "pre-ictal") %>%
  mutate(
    target_stage = factor(stage, levels = target_stage_levels),
    diff_target_minus_preictal = value - pre_value,
    metric_label_plot = factor(
      metric_bar_labels[as.character(metric)],
      levels = metric_bar_labels[metric_keys]
    )
  )

paired_summary <- diff_long %>%
  group_by(metric, metric_label, target_stage) %>%
  summarise(
    n_patients = sum(is.finite(diff_target_minus_preictal)),
    mean_diff = mean(diff_target_minus_preictal, na.rm = TRUE),
    sd_diff = sd(diff_target_minus_preictal, na.rm = TRUE),
    sem_diff = sem(diff_target_minus_preictal),
    median_diff = median(diff_target_minus_preictal, na.rm = TRUE),
    q1_diff = quantile(diff_target_minus_preictal, 0.25, na.rm = TRUE),
    q3_diff = quantile(diff_target_minus_preictal, 0.75, na.rm = TRUE),
    p_signrank = suppressWarnings(
      wilcox.test(diff_target_minus_preictal, mu = 0, exact = TRUE)$p.value
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_fdr_bh = p.adjust(p_signrank, method = "BH"),
    metric_label_plot = factor(
      metric_bar_labels[as.character(metric)],
      levels = metric_bar_labels[metric_keys]
    ),
    fdr_sig = p_fdr_bh < 0.05,
    sig_label = ifelse(fdr_sig, "*", ""),
    sig_y = mean_diff + sem_diff + 0.025
  )

scaled_patient_stage <- patient_stage
z_cols <- paste0(metric_keys, "_z")
for (i in seq_along(metric_keys)) {
  metric <- metric_keys[i]
  mu <- mean(patient_stage[[metric]], na.rm = TRUE)
  sigma <- sd(patient_stage[[metric]], na.rm = TRUE)
  scaled_patient_stage[[z_cols[i]]] <- (patient_stage[[metric]] - mu) / sigma
}

pre_scaled <- scaled_patient_stage %>%
  filter(stage == "pre-ictal") %>%
  select(patient_id, all_of(z_cols))
names(pre_scaled)[-1] <- paste0("pre_", z_cols)

state_distance <- scaled_patient_stage %>%
  filter(stage != "pre-ictal") %>%
  left_join(pre_scaled, by = "patient_id") %>%
  mutate(
    target_stage = factor(stage, levels = target_stage_levels),
    stage_num = 5 - as.numeric(target_stage)
  )

state_distance$state_distance_z <- sqrt(rowSums(
  (as.matrix(state_distance[, z_cols]) -
     as.matrix(state_distance[, paste0("pre_", z_cols)]))^2
))

state_distance_summary <- state_distance %>%
  group_by(target_stage, stage_num) %>%
  summarise(
    n_patients = sum(is.finite(state_distance_z)),
    mean_distance_z = mean(state_distance_z, na.rm = TRUE),
    sd_distance_z = sd(state_distance_z, na.rm = TRUE),
    sem_distance_z = sem(state_distance_z),
    median_distance_z = median(state_distance_z, na.rm = TRUE),
    q1_distance_z = quantile(state_distance_z, 0.25, na.rm = TRUE),
    q3_distance_z = quantile(state_distance_z, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Panel (a): patient-level repeated-measures trajectories plus mean +/- SEM.
p3a_data <- long_metrics %>%
  filter(metric == "S") %>%
  mutate(stage_num = as.numeric(stage))

p3a_summary <- stage_desc %>%
  filter(metric == "S") %>%
  mutate(stage_num = as.numeric(stage))

p3a <- ggplot(p3a_data, aes(x = stage_num, y = value)) +
  geom_line(
    aes(group = patient_id),
    colour = paper1_alpha(paper1_palette$annotation_gray, 0.25),
    linewidth = 0.24,
    show.legend = FALSE
  ) +
  geom_point(
    colour = paper1_alpha(paper1_palette$annotation_gray, 0.38),
    size = 0.78,
    show.legend = FALSE
  ) +
  geom_line(
    data = p3a_summary,
    aes(x = stage_num, y = mean, group = 1),
    inherit.aes = FALSE,
    colour = paper1_palette$primary_teal,
    linewidth = 0.68,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = p3a_summary,
    aes(x = stage_num, ymin = mean - sem, ymax = mean + sem, colour = stage),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.38,
    show.legend = FALSE
  ) +
  geom_point(
    data = p3a_summary,
    aes(x = stage_num, y = mean, colour = stage),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    size = 2.05,
    stroke = 0.64,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = stage_palette_fig3[stage_levels]) +
  scale_x_continuous(
    breaks = seq_along(stage_levels),
    labels = unname(stage_labels_short[stage_levels]),
    expand = expansion(mult = c(0.04, 0.05))
  ) +
  labs(x = NULL, y = "S") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

# Panel (a) alternative for the shallow full-width layout: patient-level
# distributions plus median/IQR and mean +/- SEM, without trajectory lines.
p3a_distribution <- ggplot(p3a_data, aes(x = stage, y = value, colour = stage)) +
  geom_boxplot(
    aes(group = stage),
    width = 0.30,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.38,
    show.legend = FALSE
  ) +
  geom_point(
    position = position_jitter(width = 0.115, height = 0, seed = 260710),
    size = 1.02,
    stroke = 0,
    alpha = 0.58,
    show.legend = FALSE
  ) +
  geom_boxplot(
    aes(group = stage),
    width = 0.30,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.38,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = p3a_summary,
    aes(x = stage, ymin = mean - sem, ymax = mean + sem, colour = stage),
    inherit.aes = FALSE,
    width = 0.09,
    linewidth = 0.42,
    show.legend = FALSE
  ) +
  geom_point(
    data = p3a_summary,
    aes(x = stage, y = mean, colour = stage),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    size = 2.05,
    stroke = 0.66,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = stage_palette_fig3[stage_levels]) +
  scale_x_discrete(
    labels = stage_labels_short,
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  labs(x = NULL, y = "S") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

# Panel (a) alternative with higher information density but no large colour
# blocks: a patient-by-stage dot heatmap plus an aligned stage summary.
p3a_patient_order <- p3a_data %>%
  filter(as.character(stage) == "pre-ictal") %>%
  arrange(value, patient_id) %>%
  pull(patient_id)

p3a_dotmap_data <- p3a_data %>%
  mutate(
    patient_rank = as.numeric(factor(patient_id, levels = p3a_patient_order)),
    stage_row = match(as.character(stage), rev(stage_levels)),
    value_label = sprintf("%.2f", value)
  )

p3a_dotmap <- ggplot(p3a_dotmap_data, aes(x = patient_rank, y = stage_row)) +
  geom_vline(
    xintercept = seq(0.5, length(p3a_patient_order) + 0.5, by = 1),
    linewidth = 0.16,
    colour = paper1_palette$grid_light
  ) +
  geom_hline(
    yintercept = seq(0.5, length(stage_levels) + 0.5, by = 1),
    linewidth = 0.16,
    colour = paper1_palette$grid_light
  ) +
  geom_point(
    aes(size = value, fill = value),
    shape = 21,
    colour = paper1_alpha(paper1_palette$annotation_gray, 0.72),
    stroke = 0.28,
    alpha = 0.82,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = value_label),
    size = FIG_HEATMAP_TEXT_PT / 2.845276,
    family = "Arial",
    colour = paper1_palette$text_dark,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(3.4, 5.5), guide = "none") +
  scale_fill_gradient(
    low = "#F2F5F4",
    high = paper1_palette$primary_teal,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq_along(p3a_patient_order),
    labels = seq_along(p3a_patient_order),
    limits = c(0.5, length(p3a_patient_order) + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(stage_levels),
    labels = unname(stage_labels_short[rev(stage_levels)]),
    limits = c(0.5, length(stage_levels) + 0.5),
    expand = c(0, 0)
  ) +
  labs(
    x = "Patients ordered by preictal S",
    y = NULL,
    title = "Individual S values"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = FIG_HEATMAP_TEXT_PT),
    plot.title = element_text(size = FIG_TEXT_PT, face = "bold", margin = margin(b = 1.2)),
    plot.margin = margin(2.0, 3.5, 2.0, 4.0)
  )

p3a_summary_dotmap <- p3a_summary %>%
  mutate(stage_row = match(as.character(stage), rev(stage_levels)))

p3a_stage_summary <- ggplot(p3a_summary_dotmap, aes(y = stage_row)) +
  geom_segment(
    aes(x = mean - sem, xend = mean + sem, yend = stage_row, colour = stage),
    linewidth = 0.48,
    show.legend = FALSE
  ) +
  geom_segment(
    aes(x = median, xend = median, y = stage_row - 0.17, yend = stage_row + 0.17),
    linewidth = 0.64,
    colour = paper1_palette$text_dark,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = mean, colour = stage),
    shape = 21,
    fill = "white",
    size = 2.02,
    stroke = 0.66,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = stage_palette_fig3[stage_levels]) +
  scale_x_continuous(
    limits = c(0, 0.9),
    breaks = c(0, 0.4, 0.8),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, length(stage_levels) + 0.5),
    expand = c(0, 0)
  ) +
  labs(x = "S", y = NULL, title = "Stage summary") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.18, colour = paper1_palette$grid_light),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.title = element_text(size = FIG_TEXT_PT, face = "bold", margin = margin(b = 1.2)),
    plot.margin = margin(2.0, 4.0, 2.0, 3.5)
  )

p3a_dotmap_combined <- (p3a_dotmap | p3a_stage_summary) +
  plot_layout(widths = c(0.80, 0.20))

# Panel (a) alternative matched to the manuscript Fig. 2(b) visual language:
# white cells containing a neutral reference line, value segment, open point,
# and direct label. Patient rows are ordered by preictal S.
p3a_cell_data <- p3a_data %>%
  mutate(
    stage_x = as.numeric(stage),
    patient_rank = as.numeric(factor(patient_id, levels = p3a_patient_order)),
    patient_y = 16 - patient_rank,
    ref_x0 = stage_x - 0.35,
    ref_x1 = stage_x + 0.10,
    point_x = ref_x0 + pmin(pmax(value, 0), 0.9) / 0.9 * (ref_x1 - ref_x0),
    label_x = stage_x + 0.29,
    value_label = sprintf("%.2f", value)
  )

p3a_cell_summary <- p3a_summary %>%
  mutate(
    stage_x = as.numeric(stage),
    summary_y = 0.65,
    ref_x0 = stage_x - 0.35,
    ref_x1 = stage_x + 0.10,
    mean_x = ref_x0 + pmin(pmax(mean, 0), 0.9) / 0.9 * (ref_x1 - ref_x0),
    low_x = ref_x0 + pmin(pmax(mean - sem, 0), 0.9) / 0.9 * (ref_x1 - ref_x0),
    high_x = ref_x0 + pmin(pmax(mean + sem, 0), 0.9) / 0.9 * (ref_x1 - ref_x0),
    label_x = stage_x + 0.29,
    value_label = sprintf("%.2f", mean)
  )

p3a_cellmatrix <- ggplot() +
  geom_tile(
    data = p3a_cell_data,
    aes(x = stage_x, y = patient_y),
    width = 0.90,
    height = 0.80,
    fill = "white",
    colour = "#D6DEDE",
    linewidth = 0.30
  ) +
  geom_tile(
    data = p3a_cell_summary,
    aes(x = stage_x, y = summary_y),
    width = 0.90,
    height = 0.82,
    fill = "white",
    colour = "#BCC8C8",
    linewidth = 0.44
  ) +
  geom_hline(
    yintercept = 1.33,
    colour = "#AEB8BA",
    linewidth = 0.38
  ) +
  geom_segment(
    data = p3a_cell_data,
    aes(x = ref_x0, xend = ref_x1, y = patient_y, yend = patient_y),
    colour = "#E2E7E7",
    linewidth = 0.20,
    lineend = "round"
  ) +
  geom_segment(
    data = p3a_cell_data,
    aes(x = ref_x0, xend = ref_x0, y = patient_y - 0.13, yend = patient_y + 0.13),
    colour = "#9BA3A8",
    linewidth = 0.28,
    lineend = "round"
  ) +
  geom_segment(
    data = p3a_cell_data,
    aes(x = ref_x0, xend = point_x, y = patient_y, yend = patient_y),
    colour = paper1_palette$primary_teal,
    linewidth = 0.66,
    lineend = "round"
  ) +
  geom_point(
    data = p3a_cell_data,
    aes(x = point_x, y = patient_y),
    shape = 21,
    fill = "white",
    colour = paper1_palette$primary_teal,
    size = 1.55,
    stroke = 0.56
  ) +
  geom_text(
    data = p3a_cell_data,
    aes(x = label_x, y = patient_y, label = value_label),
    size = FIG_HEATMAP_TEXT_PT / 2.845276,
    family = "Arial",
    colour = paper1_palette$text_dark
  ) +
  geom_segment(
    data = p3a_cell_summary,
    aes(x = ref_x0, xend = ref_x1, y = summary_y, yend = summary_y),
    colour = "#E2E7E7",
    linewidth = 0.20,
    lineend = "round"
  ) +
  geom_segment(
    data = p3a_cell_summary,
    aes(x = ref_x0, xend = ref_x0, y = summary_y - 0.13, yend = summary_y + 0.13),
    colour = "#9BA3A8",
    linewidth = 0.28,
    lineend = "round"
  ) +
  geom_segment(
    data = p3a_cell_summary,
    aes(x = low_x, xend = high_x, y = summary_y, yend = summary_y, colour = stage),
    linewidth = 0.78,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(
    data = p3a_cell_summary,
    aes(x = mean_x, y = summary_y, colour = stage),
    shape = 21,
    fill = "white",
    size = 1.82,
    stroke = 0.64,
    show.legend = FALSE
  ) +
  geom_text(
    data = p3a_cell_summary,
    aes(x = label_x, y = summary_y, label = value_label),
    size = FIG_HEATMAP_TEXT_PT / 2.845276,
    family = "Arial",
    fontface = "bold",
    colour = paper1_palette$text_dark
  ) +
  scale_colour_manual(values = stage_palette_fig3[stage_levels]) +
  scale_x_continuous(
    breaks = seq_along(stage_levels),
    labels = unname(stage_labels_short[stage_levels]),
    limits = c(0.50, length(stage_levels) + 0.50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(15:2, 0.65),
    labels = c(as.character(seq_len(14)), "Mean \u00B1 SEM"),
    limits = c(0.15, 15.45),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = NULL,
    subtitle = "Patient rows ordered by preictal S; common within-cell scale: 0\u20130.9"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = FIG_TEXT_PT),
    axis.text.y = element_text(size = FIG_HEATMAP_TEXT_PT),
    plot.subtitle = element_text(
      size = FIG_TEXT_PT,
      colour = paper1_palette$annotation_gray,
      margin = margin(b = 1.4)
    ),
    plot.margin = margin(2.0, 4.0, 2.0, 4.0)
  )

# Panel (a) moderate-density alternative: five vertical distribution strips.
# Patient observations remain visible as short horizontal ticks; stage summaries
# use the full range, IQR, median and an offset mean +/- SEM marker.
p3a_vstrip_data <- p3a_data %>%
  mutate(stage_x = as.numeric(stage))

p3a_vstrip_summary <- p3a_data %>%
  group_by(stage) %>%
  summarise(
    observed_min = min(value, na.rm = TRUE),
    observed_max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    p3a_summary %>% select(stage, mean, sem, median, q1, q3),
    by = "stage"
  ) %>%
  mutate(
    stage_x = as.numeric(stage),
    mean_x = stage_x + 0.18
  )

p3a_vstrips <- ggplot() +
  geom_segment(
    data = p3a_vstrip_data,
    aes(
      x = stage_x - 0.075,
      xend = stage_x + 0.075,
      y = value,
      yend = value
    ),
    colour = paper1_alpha(paper1_palette$annotation_gray, 0.50),
    linewidth = 0.30,
    lineend = "round"
  ) +
  geom_linerange(
    data = p3a_vstrip_summary,
    aes(x = stage_x, ymin = observed_min, ymax = observed_max),
    colour = paper1_alpha(paper1_palette$annotation_gray, 0.76),
    linewidth = 0.42
  ) +
  geom_linerange(
    data = p3a_vstrip_summary,
    aes(x = stage_x, ymin = q1, ymax = q3, colour = stage),
    linewidth = 1.24,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_segment(
    data = p3a_vstrip_summary,
    aes(
      x = stage_x - 0.13,
      xend = stage_x + 0.13,
      y = median,
      yend = median
    ),
    colour = paper1_palette$text_dark,
    linewidth = 0.62,
    lineend = "round"
  ) +
  geom_errorbar(
    data = p3a_vstrip_summary,
    aes(
      x = mean_x,
      ymin = mean - sem,
      ymax = mean + sem,
      colour = stage
    ),
    width = 0.075,
    linewidth = 0.40,
    show.legend = FALSE
  ) +
  geom_point(
    data = p3a_vstrip_summary,
    aes(x = mean_x, y = mean, colour = stage),
    shape = 21,
    fill = "white",
    size = 2.02,
    stroke = 0.66,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = stage_palette_fig3[stage_levels]) +
  scale_x_continuous(
    breaks = seq_along(stage_levels),
    labels = unname(stage_labels_short[stage_levels]),
    limits = c(0.65, length(stage_levels) + 0.40),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 0.90),
    breaks = seq(0, 0.8, by = 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "S") +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(2.0, 4.0, 2.0, 4.0)
  )

# Panel (a) balanced alternative: patient-specific paired vertical bars within
# each target-stage cluster, with compact group summaries offset to the right.
p3a_patient_offsets <- setNames(
  seq(-0.17, 0.07, length.out = length(p3a_patient_order)),
  p3a_patient_order
)

p3a_pre_values <- p3a_data %>%
  filter(as.character(stage) == "pre-ictal") %>%
  transmute(
    patient_id,
    pre_value = value,
    patient_offset = unname(p3a_patient_offsets[patient_id]),
    pair_x = 1 + patient_offset
  )

p3a_pair_data <- p3a_data %>%
  filter(as.character(stage) != "pre-ictal") %>%
  left_join(
    p3a_pre_values %>% select(patient_id, pre_value, patient_offset),
    by = "patient_id"
  ) %>%
  mutate(
    target_stage = as.character(stage),
    stage_x = as.numeric(stage),
    pair_x = stage_x + patient_offset
  )

p3a_pair_summary <- p3a_summary %>%
  mutate(
    stage_x = as.numeric(stage),
    iqr_x = stage_x + 0.17,
    mean_x = stage_x + 0.30
  )

p3a_pairbars <- ggplot() +
  geom_segment(
    data = p3a_pair_data,
    aes(x = pair_x, xend = pair_x, y = pre_value, yend = value),
    colour = paper1_alpha(paper1_palette$annotation_gray, 0.48),
    linewidth = 0.27,
    lineend = "round"
  ) +
  geom_point(
    data = p3a_pre_values,
    aes(x = pair_x, y = pre_value),
    shape = 21,
    fill = "white",
    colour = stage_palette_fig3[["pre-ictal"]],
    size = 0.96,
    stroke = 0.38,
    alpha = 0.82
  ) +
  geom_point(
    data = p3a_pair_data,
    aes(x = pair_x, y = value, colour = stage),
    shape = 21,
    fill = "white",
    size = 0.96,
    stroke = 0.40,
    alpha = 0.88,
    show.legend = FALSE
  ) +
  geom_linerange(
    data = p3a_pair_summary,
    aes(x = iqr_x, ymin = q1, ymax = q3, colour = stage),
    linewidth = 1.18,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_segment(
    data = p3a_pair_summary,
    aes(
      x = iqr_x - 0.055,
      xend = iqr_x + 0.055,
      y = median,
      yend = median
    ),
    colour = paper1_palette$text_dark,
    linewidth = 0.58,
    lineend = "round"
  ) +
  geom_errorbar(
    data = p3a_pair_summary,
    aes(
      x = mean_x,
      ymin = mean - sem,
      ymax = mean + sem,
      colour = stage
    ),
    width = 0.065,
    linewidth = 0.38,
    show.legend = FALSE
  ) +
  geom_point(
    data = p3a_pair_summary,
    aes(x = mean_x, y = mean, colour = stage),
    shape = 21,
    fill = "white",
    size = 2.00,
    stroke = 0.66,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = stage_palette_fig3[stage_levels]) +
  scale_x_continuous(
    breaks = seq_along(stage_levels),
    labels = unname(stage_labels_short[stage_levels]),
    limits = c(0.65, length(stage_levels) + 0.42),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 0.90),
    breaks = seq(0, 0.8, by = 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "S") +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(2.0, 4.0, 2.0, 4.0)
  )

# Panel (b): one vertically stacked small panel per source/sink role metric.
p3b_data <- stage_desc %>%
  filter(as.character(metric) %in% line_metric_keys) %>%
  mutate(
    metric_label_short = factor(
      metric_bar_labels[as.character(metric)],
      levels = c("D[G]", "D[H]", "D[PR]")
    ),
    stage = factor(stage, levels = stage_levels)
  )

role_extent <- max(abs(c(p3b_data$mean - p3b_data$sem, p3b_data$mean + p3b_data$sem)), na.rm = TRUE)
role_limit <- ceiling(role_extent / 0.05) * 0.05
metric_label_data <- p3b_data %>%
  distinct(metric_label_short) %>%
  mutate(
    stage = factor("pre-ictal", levels = stage_levels),
    label_y = role_limit * 0.90
  )

p3b <- ggplot(
  p3b_data,
  aes(x = stage, y = mean, colour = metric_label_short, group = metric_label_short)
) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = paper1_palette$auxiliary_gray) +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    width = 0.06,
    linewidth = 0.34,
    alpha = 0.76,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.66, alpha = 0.96, show.legend = FALSE) +
  geom_point(
    shape = 21,
    fill = "white",
    size = 1.72,
    stroke = 0.58,
    show.legend = FALSE
  ) +
  geom_text(
    data = metric_label_data,
    aes(
      x = stage,
      y = label_y,
      label = metric_label_short,
      colour = metric_label_short
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = FIG_GEOM_TEXT_SIZE,
    family = "Arial",
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_grid(rows = vars(metric_label_short), switch = "y") +
  scale_colour_manual(values = line_metric_colours_short) +
  scale_x_discrete(labels = stage_labels_short, expand = expansion(mult = c(0.04, 0.06))) +
  scale_y_continuous(limits = c(-role_limit, role_limit), breaks = seq(-role_limit, role_limit, by = 0.10)) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing.y = unit(3.0, "pt"),
    strip.text.y.left = element_blank(),
    legend.position = "none",
    plot.margin = margin(3, 4, 3, 2)
  )

make_p3b_wide_metric <- function(metric_name) {
  plot_data <- p3b_data %>% filter(as.character(metric_label_short) == metric_name)
  ggplot(plot_data, aes(x = stage, y = mean, colour = metric_label_short, group = metric_label_short)) +
    geom_hline(yintercept = 0, linewidth = 0.28, colour = paper1_palette$auxiliary_gray) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      width = 0.06,
      linewidth = 0.34,
      alpha = 0.76,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.70, alpha = 0.96, show.legend = FALSE) +
    geom_point(
      shape = 21,
      fill = "white",
      size = 1.86,
      stroke = 0.62,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = line_metric_colours_short) +
    scale_x_discrete(labels = stage_labels_short, expand = expansion(mult = c(0.03, 0.04))) +
    scale_y_continuous(
      limits = c(-role_limit, role_limit),
      breaks = c(-0.20, 0, 0.20),
      labels = c("-0.20", "0", "0.20")
    ) +
    labs(x = NULL, y = NULL, title = parse(text = metric_name)[[1]]) +
    theme(
      plot.title = element_text(
        size = FIG_TEXT_PT,
        face = "bold",
        colour = unname(line_metric_colours_short[[metric_name]]),
        margin = margin(b = 1.2)
      ),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(2.0, 4.0, 2.0, 4.0)
    )
}

p3b_wide <- make_p3b_wide_metric("D[G]") /
  make_p3b_wide_metric("D[H]") /
  make_p3b_wide_metric("D[PR]") +
  plot_layout(ncol = 1, heights = c(1, 1, 1))

p3b_3col <- (
  make_p3b_wide_metric("D[G]") |
    make_p3b_wide_metric("D[H]") |
    make_p3b_wide_metric("D[PR]")
) +
  plot_layout(nrow = 1, widths = c(1, 1, 1))

# Panel (c): paired shifts with exact Wilcoxon tests and BH correction.
p3c <- ggplot(paired_summary, aes(target_stage, mean_diff)) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = paper1_palette$auxiliary_gray) +
  geom_linerange(
    aes(ymin = 0, ymax = mean_diff, colour = target_stage),
    linewidth = 0.52,
    alpha = 0.90,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(ymin = mean_diff - sem_diff, ymax = mean_diff + sem_diff, colour = target_stage),
    width = 0.13,
    linewidth = 0.34,
    alpha = 0.82,
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = target_stage),
    shape = 21,
    fill = "white",
    size = 2.02,
    stroke = 0.66,
    show.legend = FALSE
  ) +
  geom_text(
    data = dplyr::filter(paired_summary, fdr_sig),
    aes(y = sig_y, label = sig_label),
    size = FIG_HEATMAP_TEXT_SIZE + 0.65,
    family = "Arial",
    fontface = "bold",
    colour = paper1_palette$text_dark,
    show.legend = FALSE
  ) +
  facet_wrap(~ metric_label_plot, ncol = 2, labeller = label_parsed) +
  scale_colour_manual(values = stage_outline_palette[target_stage_levels]) +
  scale_x_discrete(labels = stage_labels_short[target_stage_levels]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  labs(x = NULL, y = expression(bar(Delta*z)[list(s,k)])) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = FIG_TEXT_PT, face = "bold", lineheight = 0.78),
    panel.spacing = unit(4, "pt")
  )

make_horizontal_half_cloud <- function(data, width = 0.26) {
  pieces <- lapply(target_stage_levels, function(stage_name) {
    stage_data <- data[as.character(data$target_stage) == stage_name, , drop = FALSE]
    values <- stage_data$state_distance_z[is.finite(stage_data$state_distance_z)]
    if (length(unique(values)) < 2) return(NULL)
    y_base <- unique(stage_data$stage_num)[1]
    density_curve <- stats::density(values, from = min(values), to = max(values), n = 180)
    cloud_width <- density_curve$y / max(density_curve$y) * width
    data.frame(
      target_stage = factor(stage_name, levels = target_stage_levels),
      state_distance_z = c(density_curve$x, rev(density_curve$x)),
      cloud_y = c(rep(y_base, length(density_curve$x)), y_base + rev(cloud_width)),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(pieces)
}

p3d_cloud <- make_horizontal_half_cloud(state_distance)
set.seed(260710)
p3d_points <- state_distance %>%
  group_by(target_stage) %>%
  mutate(point_y = stage_num - 0.22 + runif(n(), -0.035, 0.035)) %>%
  ungroup()

p3d <- ggplot() +
  geom_polygon(
    data = p3d_cloud,
    aes(x = state_distance_z, y = cloud_y, group = target_stage, fill = target_stage),
    colour = NA,
    alpha = 0.14,
    show.legend = FALSE
  ) +
  geom_boxplot(
    data = state_distance,
    aes(
      x = state_distance_z,
      y = stage_num,
      group = target_stage,
      colour = target_stage
    ),
    orientation = "y",
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.36,
    show.legend = FALSE
  ) +
  geom_point(
    data = p3d_points,
    aes(x = state_distance_z, y = point_y, fill = target_stage),
    shape = 21,
    size = 1.02,
    stroke = 0.20,
    colour = paper1_alpha(paper1_palette$text_dark, 0.48),
    alpha = 0.72,
    show.legend = FALSE
  ) +
  geom_point(
    data = state_distance_summary,
    aes(x = mean_distance_z, y = stage_num),
    shape = 23,
    size = 2.04,
    stroke = 0.60,
    fill = "white",
    colour = paper1_palette$primary_teal,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = stage_palette_fig3[target_stage_levels]) +
  scale_colour_manual(values = stage_palette_fig3[target_stage_levels]) +
  scale_y_continuous(
    breaks = c(4, 3, 2, 1),
    labels = c("Early", "Mid", "Late", "Post"),
    limits = c(0.68, 4.32),
    expand = c(0, 0)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = expression(delta[list(p,s)]~"(z)"), y = NULL) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.18, colour = paper1_palette$grid_light),
    axis.ticks.y = element_blank()
  )

draw_fig3_compact <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 5,
    ncol = 2,
    widths = unit(c(0.503, 0.497), "npc"),
    heights = unit(c(0.062, 0.405, 0.038, 0.070, 0.425), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  left_x <- unit(1.8, "mm")
  top_y <- unit(1, "npc") - unit(0.8, "mm")
  draw_panel_title <- function(label, row, col) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = left_x, y = top_y, just = c("left", "top"), gp = title_gp)
    popViewport()
  }

  draw_panel_title("(a) Global DMWA symmetry", 1, 1)
  draw_panel_title("(b) Source-sink role polarization", 1, 2)
  print(p3a, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p3b, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  draw_panel_title("(c) Paired shifts from pre-ictal", 4, 1)
  draw_panel_title("(d) Mean individual direction-role distance", 4, 2)
  print(p3c, newpage = FALSE, vp = viewport(layout.pos.row = 5, layout.pos.col = 1))
  print(p3d, newpage = FALSE, vp = viewport(layout.pos.row = 5, layout.pos.col = 2))

  popViewport()
}

draw_fig3_bwide <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 6,
    ncol = 2,
    widths = unit(c(0.503, 0.497), "npc"),
    heights = unit(c(0.040, 0.250, 0.040, 0.350, 0.040, 0.280), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  left_x <- unit(1.8, "mm")
  top_y <- unit(1, "npc") - unit(0.8, "mm")
  draw_panel_title <- function(label, row, col = 1:2) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = left_x, y = top_y, just = c("left", "top"), gp = title_gp)
    popViewport()
  }

  draw_panel_title("(a) Global DMWA symmetry", 1)
  print(p3a, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))

  draw_panel_title("(b) Source-sink role polarization", 3)
  print(p3b_wide, newpage = FALSE, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))

  draw_panel_title("(c) Paired shifts from pre-ictal", 5, 1)
  draw_panel_title("(d) Mean individual direction-role distance", 5, 2)
  print(p3c, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  print(p3d, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

  popViewport()
}

draw_fig3_b3col <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 6,
    ncol = 2,
    widths = unit(c(0.503, 0.497), "npc"),
    heights = unit(c(0.040, 0.240, 0.040, 0.300, 0.040, 0.340), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  left_x <- unit(1.8, "mm")
  top_y <- unit(1, "npc") - unit(0.8, "mm")
  draw_panel_title <- function(label, row, col = 1:2) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = left_x, y = top_y, just = c("left", "top"), gp = title_gp)
    popViewport()
  }

  draw_panel_title("(a) Global DMWA symmetry", 1)
  print(p3a, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))

  draw_panel_title("(b) Source-sink role polarization", 3)
  print(p3b_3col, newpage = FALSE, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))

  draw_panel_title("(c) Paired shifts from pre-ictal", 5, 1)
  draw_panel_title("(d) Mean individual direction-role distance", 5, 2)
  print(p3c, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  print(p3d, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

  popViewport()
}

draw_fig3_b3col_adist <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 6,
    ncol = 2,
    widths = unit(c(0.503, 0.497), "npc"),
    heights = unit(c(0.040, 0.240, 0.040, 0.300, 0.040, 0.340), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  left_x <- unit(1.8, "mm")
  top_y <- unit(1, "npc") - unit(0.8, "mm")
  draw_panel_title <- function(label, row, col = 1:2) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = left_x, y = top_y, just = c("left", "top"), gp = title_gp)
    popViewport()
  }

  draw_panel_title("(a) Global DMWA symmetry", 1)
  print(p3a_distribution, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))

  draw_panel_title("(b) Source-sink role polarization", 3)
  print(p3b_3col, newpage = FALSE, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))

  draw_panel_title("(c) Paired shifts from pre-ictal", 5, 1)
  draw_panel_title("(d) Mean individual direction-role distance", 5, 2)
  print(p3c, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  print(p3d, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

  popViewport()
}

draw_fig3_b3col_adotmap <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 6,
    ncol = 2,
    widths = unit(c(0.503, 0.497), "npc"),
    heights = unit(c(0.040, 0.240, 0.040, 0.300, 0.040, 0.340), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  left_x <- unit(1.8, "mm")
  top_y <- unit(1, "npc") - unit(0.8, "mm")
  draw_panel_title <- function(label, row, col = 1:2) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = left_x, y = top_y, just = c("left", "top"), gp = title_gp)
    popViewport()
  }

  draw_panel_title("(a) Global DMWA symmetry", 1)
  print(p3a_dotmap_combined, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))

  draw_panel_title("(b) Source-sink role polarization", 3)
  print(p3b_3col, newpage = FALSE, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))

  draw_panel_title("(c) Paired shifts from pre-ictal", 5, 1)
  draw_panel_title("(d) Mean individual direction-role distance", 5, 2)
  print(p3c, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  print(p3d, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

  popViewport()
}

draw_fig3_b3col_acellmatrix <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 6,
    ncol = 2,
    widths = unit(c(0.503, 0.497), "npc"),
    heights = unit(c(0.035, 0.335, 0.035, 0.270, 0.035, 0.290), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  left_x <- unit(1.8, "mm")
  top_y <- unit(1, "npc") - unit(0.8, "mm")
  draw_panel_title <- function(label, row, col = 1:2) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = left_x, y = top_y, just = c("left", "top"), gp = title_gp)
    popViewport()
  }

  draw_panel_title("(a) Global DMWA symmetry", 1)
  print(p3a_cellmatrix, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))

  draw_panel_title("(b) Source-sink role polarization", 3)
  print(p3b_3col, newpage = FALSE, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))

  draw_panel_title("(c) Paired shifts from pre-ictal", 5, 1)
  draw_panel_title("(d) Mean individual direction-role distance", 5, 2)
  print(p3c, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  print(p3d, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

  popViewport()
}

draw_fig3_b3col_avstrips <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 6,
    ncol = 2,
    widths = unit(c(0.503, 0.497), "npc"),
    heights = unit(c(0.040, 0.240, 0.040, 0.300, 0.040, 0.340), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  left_x <- unit(1.8, "mm")
  top_y <- unit(1, "npc") - unit(0.8, "mm")
  draw_panel_title <- function(label, row, col = 1:2) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = left_x, y = top_y, just = c("left", "top"), gp = title_gp)
    popViewport()
  }

  draw_panel_title("(a) Global DMWA symmetry", 1)
  print(p3a_vstrips, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))

  draw_panel_title("(b) Source-sink role polarization", 3)
  print(p3b_3col, newpage = FALSE, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))

  draw_panel_title("(c) Paired shifts from pre-ictal", 5, 1)
  draw_panel_title("(d) Mean individual direction-role distance", 5, 2)
  print(p3c, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  print(p3d, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

  popViewport()
}

draw_fig3_b3col_apairbars <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 6,
    ncol = 2,
    widths = unit(c(0.503, 0.497), "npc"),
    heights = unit(c(0.040, 0.240, 0.040, 0.300, 0.040, 0.340), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  left_x <- unit(1.8, "mm")
  top_y <- unit(1, "npc") - unit(0.8, "mm")
  draw_panel_title <- function(label, row, col = 1:2) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, x = left_x, y = top_y, just = c("left", "top"), gp = title_gp)
    popViewport()
  }

  draw_panel_title("(a) Global DMWA symmetry", 1)
  print(p3a_pairbars, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))

  draw_panel_title("(b) Source-sink role polarization", 3)
  print(p3b_3col, newpage = FALSE, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))

  draw_panel_title("(c) Paired shifts from pre-ictal", 5, 1)
  draw_panel_title("(d) Mean individual direction-role distance", 5, 2)
  print(p3c, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  print(p3d, newpage = FALSE, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

  popViewport()
}

write.csv(p3a_data, file.path(candidate_root, "panel_a_patient_trajectories.csv"), row.names = FALSE)
write.csv(p3a_summary, file.path(candidate_root, "panel_a_stage_summary.csv"), row.names = FALSE)
write.csv(p3a_pair_data, file.path(candidate_root, "panel_a_paired_barcode.csv"), row.names = FALSE)
write.csv(p3b_data, file.path(candidate_root, "panel_b_stage_summaries.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(candidate_root, "panel_c_paired_summary.csv"), row.names = FALSE)
write.csv(state_distance, file.path(candidate_root, "panel_d_patient_distances.csv"), row.names = FALSE)
write.csv(state_distance_summary, file.path(candidate_root, "panel_d_stage_summary.csv"), row.names = FALSE)

draw_figure <- switch(
  layout_mode,
  "bwide" = draw_fig3_bwide,
  "b3col" = draw_fig3_b3col,
  "b3col_adist" = draw_fig3_b3col_adist,
  "b3col_adotmap" = draw_fig3_b3col_adotmap,
  "b3col_acellmatrix" = draw_fig3_b3col_acellmatrix,
  "b3col_avstrips" = draw_fig3_b3col_avstrips,
  "b3col_apairbars" = draw_fig3_b3col_apairbars,
  draw_fig3_compact
)
output_stem <- paste0("Fig_3_", candidate_variant)
base <- file.path(candidate_root, output_stem)
output_height_mm <- switch(
  layout_mode,
  "bwide" = 215,
  "b3col" = 198,
  "b3col_adist" = 198,
  "b3col_adotmap" = 198,
  "b3col_acellmatrix" = 225,
  "b3col_avstrips" = 198,
  "b3col_apairbars" = 198,
  148
)
manifest <- save_pub_r(draw_figure, base, width_mm = 183, height_mm = output_height_mm, dpi = 600)
write.csv(manifest, paste0(base, "_export_manifest.csv"), row.names = FALSE)

legend_draft <- paste(
  "Fig. 3. Global DMWA symmetry increases while source and sink roles polarize.",
  if (layout_mode == "b3col_apairbars") {
    "(a) Patient-level global DMWA symmetry S across seizure stages, displayed as paired vertical barcodes. Baseline points show preictal values; within each target-stage cluster, thin neutral segments connect each patient's preictal and target-stage S values and small coloured endpoints mark the target values. Offset summaries show the interquartile range, median and mean +/- SEM."
  } else if (layout_mode == "b3col_avstrips") {
    "(a) Patient-level global DMWA symmetry S across seizure stages, displayed as five vertical distribution strips. Short neutral ticks show individual patients, thin neutral lines show the observed range, thicker coloured segments show the interquartile range, dark horizontal bars mark the median, and offset open points with error bars show the mean +/- SEM."
  } else if (layout_mode == "b3col_acellmatrix") {
    "(a) Patient-level global DMWA symmetry S across seizure stages, displayed in the manuscript-style white-cell microplot matrix. Patient rows are ordered by preictal S and held fixed across columns. Within each cell, a common 0-0.9 reference scale, teal segment, open point and direct label show the raw S value; the final row shows the stage mean +/- SEM."
  } else if (layout_mode == "b3col_adotmap") {
    "(a) Patient-level global DMWA symmetry S across seizure stages. The restrained dot heatmap shows all raw patient-stage values, with patients ordered by preictal S and held in the same column across stages; circle size and low-saturation fill encode S, and labels give values to two decimal places. The aligned summary shows the mean +/- SEM and median for each stage."
  } else if (layout_mode == "b3col_adist") {
    "(a) Patient-level global DMWA symmetry S across seizure stages. Coloured points show individual patients, white boxes show the median and interquartile range, and open coloured points with error bars show the group mean +/- SEM."
  } else {
    "(a) Patient-level global DMWA symmetry S across seizure stages. Thin neutral lines and points show the five-stage trajectory of each patient; open coloured points and error bars show the group mean +/- SEM."
  },
  if (layout_mode == "bwide") {
    "(b) Stage-wise mean +/- SEM trajectories of the source-sink Gini gap, sink-minus-source entropy gap and sink-minus-source participation-ratio gap. Each metric occupies a full-width horizontal strip with a common scale and three labelled y-axis ticks."
  } else if (layout_mode %in% c("b3col", "b3col_adist", "b3col_adotmap", "b3col_acellmatrix", "b3col_avstrips", "b3col_apairbars")) {
    "(b) Stage-wise mean +/- SEM trajectories of the source-sink Gini gap, sink-minus-source entropy gap and sink-minus-source participation-ratio gap. The three metrics are arranged side-by-side so each trajectory retains full panel height, with a common scale and three labelled y-axis ticks."
  } else {
    "(b) Stage-wise mean +/- SEM trajectories of the source-sink Gini gap, sink-minus-source entropy gap and sink-minus-source participation-ratio gap, displayed in vertically aligned panels with a common scale."
  },
  "(c) Mean paired target-minus-pre-ictal shifts for S and the three source/sink role metrics. Vertical ranges start at zero, open points mark mean shifts, error bars show SEM and asterisks denote Benjamini-Hochberg-adjusted two-sided exact Wilcoxon signed-rank tests with adjusted *p* < 0.05 across the 16 stage-by-metric comparisons.",
  "(d) Patient-level displacement from each patient's own pre-ictal state in the reduced four-feature direction-role space. Light density silhouettes provide distributional context, points show individual patients, white boxes show the median and interquartile range and diamonds mark group means. This panel does not include motif features.",
  "All panels summarize n = 14 patients.",
  sep = " "
)
writeLines(legend_draft, paste0(base, "_legend_draft.txt"), useBytes = TRUE)

message("Rendered independent Fig. 3 candidate to ", base)
