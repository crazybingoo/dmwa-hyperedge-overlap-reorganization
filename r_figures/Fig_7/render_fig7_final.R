#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
  library(scales)
  library(ggtext)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
project_root <- Sys.getenv("PAPER1_PROJECT_ROOT", unset = repo_root)
stats_root <- file.path(project_root, "stats")
sim_root <- file.path(project_root, "simulation")
out_dir <- file.path(project_root, "nature_fig", "Fig_7")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fig_id <- "Fig_7"

model_levels <- c(
  "role_reciprocity_hybrid",
  "density_only",
  "source_concentration",
  "overlap_synchronization",
  "weight_strengthening",
  "reciprocity_only"
)
model_rank_levels <- rev(c(
  "source_concentration",
  "weight_strengthening",
  "density_only",
  "role_reciprocity_hybrid",
  "overlap_synchronization",
  "reciprocity_only"
))

model_labels <- c(
  "role_reciprocity_hybrid" = "Hybrid\nrole-reciprocity",
  "density_only" = "Density\nonly",
  "source_concentration" = "Source\nconcentration",
  "overlap_synchronization" = "Overlap\nsynchronization",
  "weight_strengthening" = "Weight\nstrengthening",
  "reciprocity_only" = "Reciprocity\nonly"
)
model_labels_short <- c(
  "role_reciprocity_hybrid" = "Hybrid role-reciprocity",
  "density_only" = "Density only",
  "source_concentration" = "Source concentration",
  "overlap_synchronization" = "Overlap synchronization",
  "weight_strengthening" = "Weight strengthening",
  "reciprocity_only" = "Reciprocity only"
)
model_labels_compact <- c(
  "role_reciprocity_hybrid" = "Hybrid role-recip.",
  "density_only" = "Density only",
  "source_concentration" = "Source conc.",
  "overlap_synchronization" = "Overlap sync.",
  "weight_strengthening" = "Weight strength.",
  "reciprocity_only" = "Reciprocity only"
)

model_palette <- c(
  "role_reciprocity_hybrid" = "#D2635E",
  "density_only" = "#5D83B5",
  "source_concentration" = "#2B7A78",
  "overlap_synchronization" = "#6EAD67",
  "weight_strengthening" = "#5A5F66",
  "reciprocity_only" = "#8D7BB8"
)

feature_labels <- c(
  "S" = "S",
  "source_sink_gini_gap" = "D_G",
  "reciprocal_fraction" = "f_recip",
  "dense_fraction" = "f_dense"
)

theme_paper1 <- function(base_size = 7.2, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 0.4, colour = "black"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.6),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = ggtext::element_markdown(size = base_size + 0.6, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "#5A5F66"),
      panel.grid.major = element_line(linewidth = 0.18, colour = "#EDF3F3"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}
theme_set(theme_paper1())
source(file.path(repo_root, "r_figures", "shared", "paper1_figure_style_config.R"))
panel_title_theme <- theme(
  plot.title = ggtext::element_markdown(
    size = FIG_PANEL_PT,
    face = "bold",
    hjust = 0,
    colour = paper1_palette$text_dark,
    margin = margin(b = 3)
  ),
  plot.title.position = "plot"
)

read_csv_safe <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 145, dpi = 600) {
  base <- file.path(out_dir, filename)
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white")
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in",
                res = 600, background = "white")
  print(plot)
  dev.off()

  ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in",
                 res = dpi, background = "white", compression = "lzw")
  print(plot)
  dev.off()

  data.frame(
    fig_id = fig_id,
    file = c(paste0(base, ".svg"), paste0(base, ".pdf"),
             paste0(base, ".png"), paste0(base, ".tiff")),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = c(NA, NA, dpi, dpi),
    stringsAsFactors = FALSE
  )
}

observed <- read_csv_safe(file.path(stats_root, "phase15_simulation_grid_observed_fit.csv")) %>%
  mutate(model = factor(model, levels = model_levels))

boot <- read_csv_safe(file.path(stats_root, "phase15_simulation_grid_bootstrap_summary.csv")) %>%
  mutate(model = factor(model, levels = model_levels))

rank_summary <- read_csv_safe(file.path(stats_root, "phase15_simulation_grid_ranking_stability.csv")) %>%
  mutate(model = factor(model, levels = model_levels))

rank_long <- rank_summary %>%
  select(model, top1_rate, top2_rate) %>%
  pivot_longer(c(top1_rate, top2_rate), names_to = "rank_metric", values_to = "rate") %>%
  mutate(
    rank_metric = recode(rank_metric, top1_rate = "Top-1", top2_rate = "Top-2"),
    rank_metric = factor(rank_metric, levels = c("Top-2", "Top-1"))
  )

fit <- read_csv_safe(file.path(sim_root, "competing_mechanism_fit_elbow.csv"))
empirical <- fit[1, c(
  "empirical_S",
  "empirical_source_sink_gini_gap",
  "empirical_reciprocal_fraction",
  "empirical_dense_fraction"
)]
names(empirical) <- c("S", "source_sink_gini_gap", "reciprocal_fraction", "dense_fraction")

best_observed <- observed %>%
  group_by(model) %>%
  arrange(scaled_euclidean_distance, desc(cosine_similarity), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

best_boot <- boot %>%
  group_by(model) %>%
  arrange(distance_median, desc(cosine_median), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

leading_model <- "source_concentration"
boot_leading <- boot %>% filter(as.character(model) == leading_model)
boot_context <- boot %>% filter(as.character(model) != leading_model)
best_boot_leading <- best_boot %>% filter(as.character(model) == leading_model)
best_boot_context <- best_boot %>% filter(as.character(model) != leading_model)
source_rank <- rank_summary %>% filter(as.character(model) == leading_model)
source_annotation <- data.frame(
  beta = best_boot_leading$beta,
  distance = best_boot_leading$distance_median,
  label = sprintf(
    "Source concentration\nbest β = %.1f\ntop-1 = %.3f",
    source_rank$best_beta_mode,
    source_rank$top1_rate
  )
)

a_line_labels <- c(
  "source_concentration" = "Source conc.",
  "weight_strengthening" = "Weight strength.",
  "density_only" = "Density only",
  "role_reciprocity_hybrid" = "Hybrid role-recip.",
  "overlap_synchronization" = "Overlap sync.",
  "reciprocity_only" = "Reciprocity only"
)
a_zoom_models <- c(
  "source_concentration",
  "weight_strengthening",
  "density_only",
  "role_reciprocity_hybrid"
)
a_zoom_data <- boot %>%
  filter(as.character(model) %in% a_zoom_models, beta <= 0.4)
a_zoom_context <- a_zoom_data %>%
  filter(as.character(model) != leading_model)
a_zoom_leading <- a_zoom_data %>%
  filter(as.character(model) == leading_model)

p7a_zoom <- ggplot() +
  geom_ribbon(
    data = a_zoom_context,
    aes(beta, ymin = distance_ci_low, ymax = distance_ci_high, fill = model, group = model),
    alpha = 0.075,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = a_zoom_context,
    aes(beta, distance_median, colour = model, group = model),
    linewidth = 0.46,
    alpha = 0.92,
    lineend = "round"
  ) +
  geom_point(
    data = a_zoom_context,
    aes(beta, distance_median, colour = model),
    shape = 21,
    size = 1.25,
    fill = "white",
    stroke = 0.32
  ) +
  geom_ribbon(
    data = a_zoom_leading,
    aes(beta, ymin = distance_ci_low, ymax = distance_ci_high, fill = model),
    alpha = 0.12,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = a_zoom_leading,
    aes(beta, distance_median, colour = model),
    linewidth = 0.78,
    alpha = 0.96,
    lineend = "round"
  ) +
  geom_point(
    data = a_zoom_leading,
    aes(beta, distance_median, colour = model),
    shape = 21,
    size = 1.42,
    fill = "white",
    stroke = 0.40
  ) +
  scale_colour_manual(values = model_palette, guide = "none") +
  scale_fill_manual(values = model_palette, guide = "none") +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4)) +
  scale_y_continuous(breaks = c(1.8, 2.0, 2.2, 2.4)) +
  coord_cartesian(xlim = c(0.08, 0.42), ylim = c(1.74, 2.55), clip = "on") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_paper1(base_size = FIG_HEATMAP_TEXT_PT * 0.80) +
  theme(
    plot.title = element_blank(),
    axis.text = element_text(size = FIG_HEATMAP_TEXT_PT * 0.72),
    axis.line = element_line(linewidth = 0.24, colour = paper1_palette$axis_dark),
    axis.ticks = element_line(linewidth = 0.22, colour = paper1_palette$axis_dark),
    panel.grid.major = element_line(linewidth = 0.12, colour = paper1_palette$grid_light),
    panel.border = element_rect(fill = NA, colour = paper1_palette$annotation_gray, linewidth = 0.30),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.margin = margin(2.5, 3.0, 2.5, 3.0)
  )

p7a <- ggplot() +
  geom_ribbon(
    data = boot_context,
    aes(beta, ymin = distance_ci_low, ymax = distance_ci_high, group = model, fill = model),
    alpha = 0.055, colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = boot_context,
    aes(beta, distance_median, group = model, colour = model),
    linewidth = 0.42, alpha = 0.78,
    lineend = "round"
  ) +
  geom_point(
    data = boot_context,
    aes(beta, distance_median, colour = model),
    shape = 21,
    size = 1.50,
    fill = "white",
    stroke = 0.30,
    alpha = 0.86,
    show.legend = TRUE
  ) +
  geom_ribbon(
    data = boot_leading,
    aes(beta, ymin = distance_ci_low, ymax = distance_ci_high, fill = model),
    alpha = 0.12,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = boot_leading,
    aes(beta, distance_median, colour = model),
    linewidth = 0.78,
    alpha = 0.96,
    lineend = "round"
  ) +
  geom_point(
    data = boot_leading,
    aes(beta, distance_median, colour = model),
    shape = 21,
    size = 1.70,
    fill = "white",
    stroke = 0.42
  ) +
  geom_point(
    data = best_boot_leading,
    aes(beta, distance_median),
    shape = 21,
    size = 2.35,
    fill = "white",
    colour = paper1_palette$text_dark,
    stroke = 0.42
  ) +
  geom_vline(
    xintercept = best_boot_leading$beta,
    colour = paper1_palette$annotation_gray,
    linewidth = 0.28,
    linetype = "22"
  ) +
  annotate(
    "rect",
    xmin = 0.08,
    xmax = 0.42,
    ymin = 1.74,
    ymax = 2.55,
    fill = NA,
    colour = paper1_palette$annotation_gray,
    linewidth = 0.32
  ) +
  annotate(
    "segment",
    x = 0.555,
    xend = 0.410,
    y = 9.6,
    yend = 2.55,
    colour = paper1_palette$annotation_gray,
    linewidth = 0.28,
    arrow = arrow(length = unit(3.0, "pt"), type = "closed")
  ) +
  annotate(
    "text",
    x = 0.680,
    y = 6.20,
    label = "Source concentration",
    hjust = 0,
    size = FIG_HEATMAP_TEXT_SIZE * 0.92,
    family = "Arial",
    colour = paper1_palette$text_dark
  ) +
  ggtext::geom_richtext(
    data = data.frame(
      x = 0.680,
      y = 5.45,
      label = sprintf("best <i>&beta;</i> = %.1f", source_rank$best_beta_mode)
    ),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    fill = NA,
    label.color = NA,
    label.padding = unit(c(0, 0, 0, 0), "pt"),
    label.margin = unit(c(0, 0, 0, 0), "pt"),
    size = FIG_HEATMAP_TEXT_SIZE * 0.92,
    family = "Arial",
    colour = paper1_palette$text_dark
  ) +
  annotate(
    "text",
    x = 0.680,
    y = 4.75,
    label = sprintf("top-1 = %.3f", source_rank$top1_rate),
    hjust = 0,
    size = FIG_HEATMAP_TEXT_SIZE * 0.92,
    family = "Arial",
    colour = paper1_palette$text_dark
  ) +
  annotate(
    "text",
    x = -0.120,
    y = 8.8,
    label = "D[list(m,beta)]",
    parse = TRUE,
    angle = 90,
    size = FIG_GEOM_TEXT_SIZE,
    family = "Arial",
    colour = paper1_palette$text_dark
  ) +
  scale_colour_manual(
    values = model_palette,
    breaks = names(a_line_labels),
    labels = unname(a_line_labels),
    name = NULL
  ) +
  scale_fill_manual(values = model_palette, guide = "none") +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1.0)) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1),
    breaks = c(2, 5, 10, 20, 40)
  ) +
  coord_cartesian(xlim = c(0.08, 1.02), ylim = c(1.65, 55), clip = "off") +
  labs(
    title = "(a) Distance to empirical target",
    x = "Mechanism strength <i>&beta;</i>",
    y = NULL
  ) +
  guides(
    colour = guide_legend(
      ncol = 2,
      byrow = TRUE,
      override.aes = list(
        linewidth = 0.82,
        alpha = 1,
        size = 2.05,
        shape = 21,
        fill = "white",
        stroke = 0.42
      )
    )
  ) +
  theme(
    legend.position = c(0.020, 1.035),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = paper1_alpha("#FFFFFF", 0.82), colour = NA),
    legend.key.width = unit(14.0, "pt"),
    legend.key.height = unit(6.4, "pt"),
    legend.spacing.x = unit(1.4, "pt"),
    legend.spacing.y = unit(0.1, "pt"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(1, 1, 1, 1),
    legend.text = element_text(size = FIG_HEATMAP_TEXT_PT * 0.92, colour = paper1_palette$text_dark),
    panel.grid.major.x = element_line(linewidth = 0.16, colour = paper1_palette$grid_light),
    panel.grid.major.y = element_line(linewidth = 0.16, colour = paper1_palette$grid_light),
    axis.text.y = element_text(margin = margin(r = 0.2)),
    axis.title.x = ggtext::element_markdown(),
    plot.margin = margin(5, 1, 5, 31)
  ) +
  panel_title_theme

p7a <- p7a + inset_element(
  p7a_zoom,
  left = 0.462,
  bottom = 0.455,
  right = 0.907,
  top = 0.765,
  align_to = "panel",
  on_top = TRUE,
  clip = TRUE
)

rank_plot <- rank_long %>%
  mutate(
    model = factor(as.character(model), levels = model_rank_levels),
    model_y = as.numeric(model),
    rank_offset = ifelse(rank_metric == "Top-1", 0.18, -0.18),
    y_plot = model_y + rank_offset,
    label = sprintf("%.3f", rate)
  )
rank_label_plot <- rank_plot %>%
  filter(
    (as.character(model) == leading_model & rank_metric == "Top-1") |
      (as.character(model) == "weight_strengthening" & rank_metric == "Top-2") |
      (as.character(model) == "density_only" & rank_metric == "Top-2")
  )

p7b <- ggplot(rank_plot, aes(rate, y_plot, colour = rank_metric)) +
  geom_segment(aes(x = 0, xend = rate, yend = y_plot),
               linewidth = 0.62, lineend = "round", alpha = 0.92) +
  geom_point(shape = 21, size = 2.05, fill = "white", stroke = 0.38) +
  geom_text(
    data = rank_label_plot,
    aes(x = pmin(rate + 0.035, 1.03), label = label),
    hjust = 0,
    size = FIG_HEATMAP_TEXT_SIZE * 0.90,
    family = "Arial",
    colour = paper1_palette$annotation_gray,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Top-2" = paper1_palette$auxiliary_gray,
                                 "Top-1" = paper1_palette$primary_teal),
                      name = NULL) +
  scale_y_continuous(
    breaks = seq_along(model_rank_levels),
    labels = model_labels_compact[model_rank_levels],
    limits = c(0.45, length(model_rank_levels) + 0.55),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(limits = c(0, 1.10), breaks = c(0, 0.5, 1.0),
                     labels = number_format(accuracy = 0.1)) +
  labs(
    title = "(b) Ranking stability",
    x = expression(pi[m]^{(1)}*","~pi[m]^{("<=2")}),
    y = NULL
  ) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.8, size = 2.1))) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.16, colour = paper1_palette$grid_light),
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = unit(11, "pt"),
    legend.spacing.x = unit(4, "pt"),
    legend.box.margin = margin(-4, 0, -7, 0),
    axis.text.y = element_text(size = FIG_TEXT_PT - 1.75, lineheight = 0.96),
    plot.margin = margin(5, 4, 4, 7)
  ) +
  panel_title_theme

sensitivity <- read_csv_safe(file.path(stats_root, "results_logic_simulation_patient_target_sensitivity.csv"))

original_best <- observed %>%
  group_by(model) %>%
  slice_min(scaled_euclidean_distance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(model = as.character(model), variant = "Original scaled",
            beta = beta, score = scaled_euclidean_distance)

sensitivity_long <- bind_rows(
  original_best,
  sensitivity %>% transmute(model, variant = "Patient unscaled",
                            beta = best_beta_unscaled, score = min_unscaled_distance),
  sensitivity %>% transmute(model, variant = "Patient no dense",
                            beta = best_beta_no_dense, score = min_unscaled_distance_no_dense),
  sensitivity %>% transmute(model, variant = "Patient cosine",
                            beta = best_beta_cosine, score = -max_cosine_similarity)
) %>%
  group_by(variant) %>%
  mutate(rank = rank(score, ties.method = "min")) %>%
  ungroup() %>%
  mutate(
    model_label = factor(model_labels_compact[model], levels = model_labels_compact[model_rank_levels]),
    variant = factor(variant, levels = c("Original scaled", "Patient unscaled", "Patient no dense", "Patient cosine")),
    rank_label = as.character(rank),
    beta_label = sprintf("<i>&beta;</i> = %.1f", beta)
  )

p7c_data <- sensitivity_long %>%
  mutate(
    variant_index = as.numeric(variant),
    model_short = unname(model_labels_compact[model]),
    emphasis = model %in% c("source_concentration", "role_reciprocity_hybrid")
  )

p7c_context <- p7c_data %>% filter(!emphasis)
p7c_emphasis <- p7c_data %>% filter(emphasis)
p7c_end_context <- p7c_context %>% filter(variant_index == max(variant_index))
p7c_end_emphasis <- p7c_emphasis %>% filter(variant_index == max(variant_index))

p7c <- ggplot(p7c_data, aes(variant_index, rank, group = model)) +
  geom_hline(
    yintercept = 1:6,
    linewidth = 0.16,
    colour = paper1_palette$grid_light
  ) +
  geom_line(
    data = p7c_context,
    aes(colour = model),
    linewidth = 0.38,
    alpha = 0.46,
    show.legend = FALSE
  ) +
  geom_point(
    data = p7c_context,
    aes(colour = model),
    shape = 21,
    fill = "white",
    size = 1.50,
    stroke = 0.38,
    alpha = 0.72,
    show.legend = FALSE
  ) +
  geom_line(
    data = p7c_emphasis,
    aes(colour = model),
    linewidth = 0.74,
    alpha = 0.98,
    show.legend = FALSE
  ) +
  geom_point(
    data = p7c_emphasis,
    aes(colour = model),
    shape = 21,
    fill = "white",
    size = 2.20,
    stroke = 0.62,
    show.legend = FALSE
  ) +
  geom_text(
    data = p7c_end_context,
    aes(x = 4.18, label = model_short, colour = model),
    hjust = 0,
    size = FIG_HEATMAP_TEXT_SIZE * 0.78,
    family = "Arial",
    alpha = 0.74,
    show.legend = FALSE
  ) +
  geom_text(
    data = p7c_end_emphasis,
    aes(x = 4.18, label = model_short, colour = model),
    hjust = 0,
    size = FIG_HEATMAP_TEXT_SIZE * 0.84,
    family = "Arial",
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_colour_manual(values = model_palette, guide = "none") +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("Original\nscaled", "Patient\nunscaled",
               "Patient\nno dense", "Patient\ncosine"),
    limits = c(0.92, 5.20),
    expand = expansion(mult = 0)
  ) +
  scale_y_reverse(
    breaks = 1:6,
    limits = c(6.25, 0.75),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(title = "(c) Rank trajectories", x = NULL, y = expression(r[m])) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(linewidth = 0.30, colour = paper1_palette$axis_dark),
    axis.line.y = element_line(linewidth = 0.30, colour = paper1_palette$axis_dark),
    axis.ticks.x = element_line(linewidth = 0.26, colour = paper1_palette$axis_dark),
    axis.ticks.y = element_line(linewidth = 0.26, colour = paper1_palette$axis_dark),
    axis.text.x = element_text(size = FIG_TEXT_PT - 1.80, lineheight = 0.92),
    axis.text.y = element_text(size = FIG_TEXT_PT - 1.90),
    axis.title.y = element_text(size = FIG_TEXT_PT - 1.50, margin = margin(r = 1.5)),
    plot.margin = margin(5, 3, 4, 0)
  ) +
  panel_title_theme

component_ratios <- best_observed %>%
  select(model, beta, starts_with("simulated_")) %>%
  pivot_longer(starts_with("simulated_"), names_to = "feature", values_to = "simulated") %>%
  mutate(
    feature = sub("^simulated_", "", feature),
    empirical = as.numeric(empirical[feature]),
    ratio = simulated / empirical,
    model = factor(as.character(model), levels = model_rank_levels),
    feature = factor(feature, levels = names(feature_labels)),
    model_y = as.numeric(model),
    ratio_label = sprintf("%.2f", ratio)
  ) %>%
  group_by(feature) %>%
  mutate(
    feature_span = max(ratio, na.rm = TRUE) - min(ratio, na.rm = TRUE),
    label_pad = pmax(feature_span * 0.06, 0.018),
    ratio_label_x = ifelse(ratio >= 0, ratio + label_pad, ratio - label_pad),
    ratio_label_hjust = ifelse(ratio >= 0, 0, 1)
  ) %>%
  ungroup()

d_row_bands <- data.frame(
  feature = rep(factor(names(feature_labels), levels = names(feature_labels)),
                each = length(model_rank_levels)),
  row_id = seq_along(model_rank_levels)
) %>%
  mutate(
    ymin = row_id - 0.5,
    ymax = row_id + 0.5
  ) %>%
  filter(row_id %% 2 == 1)

feature_lollipop_palette <- c(
  "S" = paper1_palette$pre,
  "source_sink_gini_gap" = paper1_palette$primary_teal,
  "reciprocal_fraction" = paper1_palette$accent_red,
  "dense_fraction" = paper1_palette$post
)

feature_labels_d <- c(
  "S" = "S",
  "source_sink_gini_gap" = "D<sub>G</sub>",
  "reciprocal_fraction" = "f<sub>recip</sub>",
  "dense_fraction" = "f<sub>dense</sub>"
)

d_ratio_reference <- component_ratios %>%
  group_by(feature) %>%
  summarise(
    ratio_min = min(ratio, na.rm = TRUE),
    ratio_max = max(ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(1 >= ratio_min & 1 <= ratio_max)

p7d <- ggplot(component_ratios, aes(ratio, model_y)) +
  geom_rect(
    data = d_row_bands,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = paper1_alpha(paper1_palette$grid_light, 0.44),
    colour = NA
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.26,
    colour = paper1_palette$auxiliary_gray
  ) +
  geom_vline(
    data = d_ratio_reference,
    aes(xintercept = 1),
    inherit.aes = FALSE,
    linewidth = 0.24,
    colour = paper1_palette$annotation_gray,
    linetype = "22"
  ) +
  geom_segment(
    aes(x = 0, xend = ratio, yend = model_y, colour = feature),
    linewidth = 0.48,
    alpha = 0.82,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = feature),
    shape = 21,
    size = 2.05,
    fill = "white",
    stroke = 0.42,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = ratio_label_x, label = ratio_label, hjust = ratio_label_hjust),
    size = FIG_HEATMAP_TEXT_SIZE * 0.90,
    family = "Arial",
    colour = paper1_palette$annotation_gray
  ) +
  scale_colour_manual(values = feature_lollipop_palette, guide = "none") +
  scale_x_continuous(
    breaks = scales::breaks_extended(n = 3),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.16, 0.28))
  ) +
  scale_y_continuous(
    breaks = seq_along(model_rank_levels),
    labels = model_labels_compact[model_rank_levels],
    expand = expansion(mult = c(0, 0), add = c(0.32, 0.32))
  ) +
  facet_grid(. ~ feature, scales = "free_x", labeller = labeller(feature = feature_labels_d)) +
  coord_cartesian(ylim = c(0.5, length(model_rank_levels) + 0.5), clip = "off") +
  labs(
    title = "(d) Component reproduction",
    x = expression(R[list(m,beta,k)] == "simulated / empirical"),
    y = NULL
  ) +
  theme(
    strip.text = ggtext::element_markdown(
      size = FIG_TEXT_PT - 1.0,
      face = "bold",
      colour = paper1_palette$text_dark,
      lineheight = 0.92
    ),
    axis.text.x = element_text(size = FIG_TEXT_PT - 0.45),
    axis.text.y = element_text(size = FIG_TEXT_PT - 1.65, lineheight = 0.96),
    axis.title.x = element_text(margin = margin(t = 1.5)),
    panel.grid = element_blank(),
    panel.spacing.x = unit(6, "pt"),
    legend.position = "none",
    plot.margin = margin(7, 3, 4, 4)
  ) +
  panel_title_theme

fig7_design <- "
AAB
AAC
DDD
"

fig7 <- free(p7a, type = "panel", side = "l") + p7b + free(p7c, type = "panel", side = "l") + p7d +
  plot_layout(
    design = fig7_design,
    widths = c(0.84, 0.50, 1.56),
    heights = c(1.05, 1.05, 0.98),
    guides = "keep"
  ) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

manifest <- save_pub_r(fig7, fig_id, 183, 158)
manifest$exists <- file.exists(manifest$file)
manifest$size_bytes <- file.info(manifest$file)$size
write.csv(manifest, file.path(out_dir, paste0(fig_id, "_export_manifest.csv")), row.names = FALSE)

source_data <- bind_rows(
  boot %>%
    transmute(panel = "a_grid_distance", model = as.character(model), beta,
              value = distance_median, ci_low = distance_ci_low,
              ci_high = distance_ci_high, feature = NA_character_),
  rank_long %>%
    transmute(panel = "b_ranking_rate_bar", model = as.character(model), beta = NA_real_,
              value = rate, ci_low = NA_real_, ci_high = NA_real_,
              feature = as.character(rank_metric)),
  sensitivity_long %>%
    transmute(panel = "c_rank_trajectories", model = as.character(model),
              beta, value = rank, ci_low = NA_real_, ci_high = NA_real_,
              feature = as.character(variant)),
  component_ratios %>%
    transmute(panel = "d_component_reproduction", model = as.character(model),
              beta, value = ratio, ci_low = NA_real_, ci_high = NA_real_,
              feature = as.character(feature))
)
write.csv(source_data, file.path(out_dir, paste0(fig_id, "_source_data.csv")), row.names = FALSE)

beta_char <- intToUtf8(0x03B2)
caption <- paste(
  "Fig. 7. Full simulation-grid evaluation of competing structural mechanisms.",
  paste0("a, Bootstrap median scaled Euclidean distance between each simulated four-feature effect vector and the empirical ictal effect vector across the ", beta_char, " grid. Coloured lines and points denote individual mechanisms, shaded bands show 95% bootstrap intervals, and the inset magnifies the low-distance range from ", beta_char, " = 0.1 to 0.4."),
  paste0("b, Frequency with which each model was ranked first or second after bootstrap resampling and ", beta_char, " selection, shown as horizontal selection-rate intervals."),
  paste0("c, Rank trajectories after ", beta_char, " selection under the original scaled seizure-level target and three patient-level target/scaling variants. Lines connect each model across variants; source concentration and the hybrid role-reciprocity model are emphasized because they exchange the leading positions as the target definition and distance measure change."),
  paste0("d, Ratio of simulated to empirical feature effects at each model's observed best ", beta_char, "."),
  "Negative ratios indicate simulated effects with the opposite sign to the empirical component.",
  sep = " "
)
writeLines(enc2utf8(caption), file.path(out_dir, paste0(fig_id, "_legend_draft.txt")), useBytes = TRUE)

message("Updated ", fig_id, " exported to: ", out_dir)
print(manifest)





