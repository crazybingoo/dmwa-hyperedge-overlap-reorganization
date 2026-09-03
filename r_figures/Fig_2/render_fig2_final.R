#!/usr/bin/env Rscript
# Variant of render_fig2_subfigure_titles_bold_ac_trimmed_R.R that
# changes the line+marker colour of subplot (a) to blue, and the
# scatter points of subplot (c) to the same blue at the same 0.42
# alpha. Keeps the rest of the figure identical. Does not overwrite
# the original Fig_2_subfigure_titles_bold_ac_trimmed.* outputs.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(svglite)
  library(ragg)
  library(grid)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
project_root <- Sys.getenv("PAPER1_PROJECT_ROOT", unset = repo_root)
source_root <- file.path(project_root, "figures_data", "Fig_2")
out_root <- file.path(project_root, "nature_fig", "Fig_2")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

stage_levels <- c("pre-ictal", "early", "mid", "late", "post-ictal")
contrast_levels <- c("early", "mid", "late", "post-ictal")

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
      plot.title = element_text(size = base_size + 0.6, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#5A5F66"),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#EDF3F3"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

theme_set(theme_paper1())
source(file.path(repo_root, "r_figures", "shared", "paper1_figure_style_config.R"))
theme_set(theme_paper1(base_size = FIG_TEXT_PT))

wrap_text <- function(x, width = 18) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

macro_levels <- c("NIZ", "PZ", "PZ+NIZ", "SOZ+NIZ / SOZ+PZ+NIZ", "SOZ / SOZ+PZ")
macro_levels_a <- wrap_text(macro_levels, 18)
macro_levels_b <- wrap_text(macro_levels, 22)

stage_factor <- function(x) {
  factor(x, levels = stage_levels)
}

contrast_stage_factor <- function(x) {
  factor(x, levels = contrast_levels)
}

stage_labels_short <- c(
  "pre-ictal" = "Pre",
  "early" = "Early",
  "mid" = "Mid",
  "late" = "Late",
  "post-ictal" = "Post"
)

read_csv_safe <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

save_pub_grid <- function(draw_fun, base, width_mm = 183, height_mm = 132, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white")
  draw_fun()
  dev.off()

  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = "Arial")
  draw_fun()
  dev.off()

  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = 600, background = "white")
  draw_fun()
  dev.off()

  ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = dpi, background = "white", compression = "lzw")
  draw_fun()
  dev.off()
}

summary_path <- file.path(source_root, "fig2_macro_state_summary_elbow.csv")
contrast_path <- file.path(source_root, "fig2_macro_state_paired_contrasts_elbow.csv")
patient_path <- file.path(source_root, "fig2_macro_state_patient_occupancy_elbow.csv")

fig2_sum <- read_csv_safe(summary_path) %>%
  mutate(
    stage = stage_factor(stage),
    macro_state = factor(macro_state, levels = macro_levels),
    macro_label = factor(wrap_text(as.character(macro_state), 18), levels = macro_levels_a)
  )

fig2_con <- read_csv_safe(contrast_path) %>%
  mutate(
    target_stage = contrast_stage_factor(target_stage),
    macro_state = factor(macro_state, levels = macro_levels),
    macro_label = factor(wrap_text(as.character(macro_state), 22), levels = rev(macro_levels_b)),
    fdr_sig = signrank_p_fdr_bh < 0.05,
    fdr_label = ifelse(signrank_p_fdr_bh < 0.05, "*", ""),
    diff_label = ifelse(
      fdr_sig,
      paste0(sprintf("%+.1f", mean_diff_target_minus_preictal), "*"),
      sprintf("%+.1f", mean_diff_target_minus_preictal)
    ),
    abs_diff = abs(mean_diff_target_minus_preictal),
    max_abs_diff = max(abs(mean_diff_target_minus_preictal), na.rm = TRUE),
    stage_x = as.numeric(target_stage),
    macro_y = as.numeric(macro_label),
    diff_dir = ifelse(mean_diff_target_minus_preictal >= 0, "Increase", "Decrease"),
    zero_x = stage_x - 0.10,
    delta_x = scales::rescale(abs_diff, to = c(0.040, 0.255), from = c(0, max(max_abs_diff, na.rm = TRUE))),
    seg_x = ifelse(mean_diff_target_minus_preictal >= 0, zero_x, zero_x - delta_x),
    seg_xend = ifelse(mean_diff_target_minus_preictal >= 0, zero_x + delta_x, zero_x),
    point_x = ifelse(mean_diff_target_minus_preictal >= 0, seg_xend, seg_x),
    value_x = stage_x + 0.29
  )

safe_wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(unique(x)) < 2) {
    return(NA_real_)
  }
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

patient_case_occ <- read_csv_safe(patient_path) %>%
  mutate(
    stage = stage_factor(stage),
    macro_state = factor(macro_state, levels = macro_levels)
  )

patient_occ <- patient_case_occ %>%
  group_by(patient_id, stage, macro_state_index, macro_state) %>%
  summarise(
    occupancy_percent = mean(occupancy_percent, na.rm = TRUE),
    n_cases = n_distinct(case_id),
    .groups = "drop"
  )

patient_occ_summary <- patient_occ %>%
  group_by(stage, macro_state_index, macro_state) %>%
  summarise(
    mean_occupancy_percent = mean(occupancy_percent, na.rm = TRUE),
    sem = sd(occupancy_percent, na.rm = TRUE) / sqrt(sum(is.finite(occupancy_percent))),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  mutate(
    macro_state = factor(macro_state, levels = macro_levels),
    macro_label = factor(wrap_text(as.character(macro_state), 18), levels = macro_levels_a)
  )

patient_pre <- patient_occ %>%
  filter(stage == "pre-ictal") %>%
  select(patient_id, macro_state_index, macro_state, preictal_occupancy_percent = occupancy_percent)

patient_target <- patient_occ %>%
  filter(as.character(stage) %in% contrast_levels) %>%
  transmute(
    patient_id,
    target_stage = as.character(stage),
    macro_state_index,
    macro_state,
    target_occupancy_percent = occupancy_percent
  )

patient_diff_long <- patient_target %>%
  inner_join(patient_pre, by = c("patient_id", "macro_state_index", "macro_state")) %>%
  mutate(
    diff_target_minus_preictal = target_occupancy_percent - preictal_occupancy_percent,
    target_stage = contrast_stage_factor(target_stage)
  )

fig2_patient_con <- patient_diff_long %>%
  group_by(target_stage, macro_state_index, macro_state) %>%
  summarise(
    n_patients = n_distinct(patient_id),
    preictal_mean = mean(preictal_occupancy_percent, na.rm = TRUE),
    target_mean = mean(target_occupancy_percent, na.rm = TRUE),
    mean_diff_target_minus_preictal = mean(diff_target_minus_preictal, na.rm = TRUE),
    median_diff_target_minus_preictal = median(diff_target_minus_preictal, na.rm = TRUE),
    positive_fraction = mean(diff_target_minus_preictal > 0, na.rm = TRUE),
    signrank_p = safe_wilcox_p(diff_target_minus_preictal),
    .groups = "drop"
  ) %>%
  mutate(
    signrank_p_fdr_bh = p.adjust(signrank_p, method = "BH"),
    fdr_sig = !is.na(signrank_p_fdr_bh) & signrank_p_fdr_bh < 0.05,
    macro_label = factor(wrap_text(as.character(macro_state), 22), levels = rev(macro_levels_b)),
    abs_diff = abs(mean_diff_target_minus_preictal),
    max_abs_diff = max(abs(mean_diff_target_minus_preictal), na.rm = TRUE),
    stage_x = as.numeric(target_stage),
    macro_y = as.numeric(macro_label),
    diff_dir = ifelse(mean_diff_target_minus_preictal >= 0, "Increase", "Decrease"),
    zero_x = stage_x - 0.10,
    delta_x = scales::rescale(abs_diff, to = c(0.040, 0.255), from = c(0, max(max_abs_diff, na.rm = TRUE))),
    seg_x = ifelse(mean_diff_target_minus_preictal >= 0, zero_x, zero_x - delta_x),
    seg_xend = ifelse(mean_diff_target_minus_preictal >= 0, zero_x + delta_x, zero_x),
    point_x = ifelse(mean_diff_target_minus_preictal >= 0, seg_xend, seg_x),
    value_x = stage_x + 0.29,
    diff_label = sprintf("%+.1f", mean_diff_target_minus_preictal)
  )

patient_ictal_diff <- patient_diff_long %>%
  filter(as.character(target_stage) %in% c("early", "mid", "late")) %>%
  group_by(patient_id, macro_state_index, macro_state) %>%
  summarise(
    ictal_mean_diff = mean(diff_target_minus_preictal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    macro_label = factor(wrap_text(as.character(macro_state), 22), levels = rev(macro_levels_b)),
    macro_y = as.numeric(macro_label)
  ) %>%
  group_by(macro_label) %>%
  arrange(patient_id, .by_group = TRUE) %>%
  mutate(point_y = macro_y + seq(-0.16, 0.16, length.out = n())) %>%
  ungroup()

patient_ictal_summary <- patient_ictal_diff %>%
  group_by(macro_label, macro_y) %>%
  summarise(
    mean_ictal_diff = mean(ictal_mean_diff, na.rm = TRUE),
    .groups = "drop"
  )

case_patient_agreement <- fig2_con %>%
  transmute(
    target_stage,
    macro_state_index,
    macro_state,
    case_mean_diff = mean_diff_target_minus_preictal
  ) %>%
  inner_join(
    fig2_patient_con %>%
      transmute(
        target_stage,
        macro_state_index,
        macro_state,
        patient_mean_diff = mean_diff_target_minus_preictal,
        patient_fdr_sig = fdr_sig
      ),
    by = c("target_stage", "macro_state_index", "macro_state")
  )

agreement_lim <- ceiling(max(abs(c(case_patient_agreement$case_mean_diff, case_patient_agreement$patient_mean_diff)), na.rm = TRUE) / 5) * 5
agreement_lim <- max(agreement_lim, 15)
agreement_r <- cor(case_patient_agreement$case_mean_diff, case_patient_agreement$patient_mean_diff, use = "complete.obs")

# Blue colour override for subplots (a) and (c) scatter
line_blue <- "#1F77B4"
line_blue_alpha42 <- paper1_alpha(line_blue, 0.42)

p2a <- ggplot(fig2_sum, aes(stage, mean_occupancy_percent, group = 1)) +
  geom_line(colour = line_blue, linewidth = 0.68, alpha = 0.94) +
  geom_errorbar(
    aes(ymin = mean_occupancy_percent - sem, ymax = mean_occupancy_percent + sem),
    width = 0.06, linewidth = 0.34, colour = line_blue, alpha = 0.62
  ) +
  geom_point(shape = 21, size = 2.10, fill = "white",
             colour = line_blue, stroke = 0.58) +
  facet_wrap(~ macro_label, nrow = 1) +
  scale_x_discrete(labels = stage_labels_short) +
  labs(y = expression(Omega[g]~"(%)"), x = NULL) +
  theme(axis.text.x = element_text(angle = 32, hjust = 1))

p2b <- ggplot(
  fig2_con,
  aes(x = stage_x, y = macro_y)
) +
  geom_tile(width = 0.88, height = 0.82, fill = "white",
            colour = "#D6DEDE", linewidth = 0.34) +
  geom_segment(
    aes(x = zero_x - 0.27, xend = zero_x + 0.27, y = macro_y, yend = macro_y),
    colour = "#E2E7E7", linewidth = 0.20, lineend = "round"
  ) +
  geom_segment(
    aes(x = zero_x, xend = zero_x, y = macro_y - 0.13, yend = macro_y + 0.13),
    colour = "#9BA3A8", linewidth = 0.30, lineend = "round"
  ) +
  geom_segment(
    aes(x = seg_x, xend = seg_xend, y = macro_y, yend = macro_y, colour = diff_dir),
    linewidth = 0.78, lineend = "round", show.legend = FALSE
  ) +
  geom_point(
    aes(x = point_x, y = macro_y, colour = diff_dir),
    shape = 21, fill = "white", size = 2.15, stroke = 0.68, show.legend = FALSE
  ) +
  geom_text(
    aes(x = value_x, label = diff_label),
    size = FIG_HEATMAP_TEXT_SIZE, family = "Arial", colour = paper1_palette$text_dark
  ) +
  scale_colour_manual(
    values = c("Decrease" = "#4C72B0", "Increase" = "#C44E52"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq_along(contrast_levels),
    labels = stage_labels_short[contrast_levels],
    limits = c(0.43, length(contrast_levels) + 0.62),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(fig2_con$macro_label)),
    labels = levels(fig2_con$macro_label),
    limits = c(0.5, length(levels(fig2_con$macro_label)) + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    y = NULL,
    x = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title.position = "plot"
  )

p2c <- ggplot(patient_ictal_diff, aes(x = ictal_mean_diff, y = macro_y, group = macro_y)) +
  geom_vline(xintercept = 0, colour = "#AEB6BA", linewidth = 0.32, linetype = "dashed") +
  geom_boxplot(
    orientation = "y", width = 0.40, fill = "white",
    colour = paper1_palette$axis_dark, linewidth = 0.36,
    outlier.shape = NA
  ) +
  geom_point(
    aes(x = ictal_mean_diff, y = point_y),
    shape = 16, size = 0.95,
    colour = line_blue_alpha42
  ) +
  geom_point(
    data = patient_ictal_summary,
    aes(x = mean_ictal_diff, y = macro_y),
    inherit.aes = FALSE,
    shape = 23, size = 2.10,
    fill = "white", colour = paper1_palette$accent_red, stroke = 0.52
  ) +
  scale_x_continuous(
    breaks = c(-50, -25, 0, 25, 50),
    limits = c(-55, 55),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(patient_ictal_diff$macro_label)),
    labels = levels(patient_ictal_diff$macro_label),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    x = expression("Ictal mean"~Delta*Omega[g]~"(pp)"),
    y = NULL
  ) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.16, colour = paper1_palette$grid_light),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = FIG_TEXT_PT - 0.3),
    plot.margin = margin(4, 4, 4, 4)
  )

p2d <- ggplot(case_patient_agreement, aes(x = case_mean_diff, y = patient_mean_diff)) +
  geom_hline(yintercept = 0, colour = "#D8DEDE", linewidth = 0.24) +
  geom_vline(xintercept = 0, colour = "#D8DEDE", linewidth = 0.24) +
  geom_abline(slope = 1, intercept = 0, colour = "#8A9296", linewidth = 0.36, linetype = "dashed") +
  geom_point(
    shape = 21, size = 2.05, stroke = 0.58,
    colour = paper1_palette$primary_teal,
    fill = "white"
  ) +
  annotate(
    "text",
    x = -agreement_lim * 0.94,
    y = agreement_lim * 0.88,
    label = sprintf("20 contrasts\n14 patients\nr = %.2f", agreement_r),
    hjust = 0, vjust = 1,
    family = "Arial", size = FIG_HEATMAP_TEXT_SIZE,
    colour = paper1_palette$annotation_gray,
    lineheight = 0.95
  ) +
  coord_equal(xlim = c(-agreement_lim, agreement_lim), ylim = c(-agreement_lim, agreement_lim), expand = FALSE) +
  labs(
    title = "(d) Case-patient agreement",
    x = expression("Seizure-level"~Delta*Omega[g]~"(pp)"),
    y = expression("Patient-level"~Delta*Omega[g]~"(pp)")
  ) +
  theme(
    plot.title = element_text(size = FIG_PANEL_PT, face = "bold", hjust = 0),
    panel.grid.major = element_blank(),
    axis.title.y = element_text(margin = margin(r = 1.5)),
    plot.margin = margin(4, 4, 4, 6)
  )

draw_fig2 <- function() {
  grid.newpage()
  layout <- grid.layout(
    nrow = 8,
    ncol = 1,
    heights = unit(c(0.052, 0.300, 0.024, 0.098, 0.265, 0.032, 0.052, 0.177), "npc")
  )
  pushViewport(viewport(layout = layout))

  title_gp <- gpar(
    fontfamily = "Arial",
    fontface = "bold",
    fontsize = FIG_PANEL_PT,
    col = paper1_palette$text_dark
  )
  subtitle_gp <- gpar(
    fontfamily = "Arial",
    fontsize = FIG_TEXT_PT,
    col = paper1_palette$annotation_gray,
    lineheight = 1.05
  )
  left_x <- unit(1.8, "mm")

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(
    "(a) Macro-state occupancy",
    x = left_x, y = unit(1, "npc") - unit(0.8, "mm"),
    just = c("left", "top"), gp = title_gp
  )
  popViewport()

  print(p2a, newpage = FALSE, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
  grid.text(
    "(b) Seizure-level mean change from pre-ictal",
    x = left_x, y = unit(1, "npc") - unit(0.6, "mm"),
    just = c("left", "top"), gp = title_gp
  )
  grid.text(
    expression("Values are mean paired "*Delta*Omega[g]*" (pp);"),
    x = left_x, y = unit(1, "npc") - unit(5.0, "mm"),
    just = c("left", "top"), gp = subtitle_gp
  )
  grid.text(
    "* denotes paired Wilcoxon, Benjamini-Hochberg FDR < 0.05",
    x = left_x, y = unit(1, "npc") - unit(9.2, "mm"),
    just = c("left", "top"), gp = subtitle_gp
  )
  popViewport()

  print(p2b, newpage = FALSE, vp = viewport(layout.pos.row = 5, layout.pos.col = 1))

  pushViewport(viewport(layout.pos.row = 7, layout.pos.col = 1))
  grid.text(
    "(c) Patient-level ictal shift",
    x = left_x, y = unit(1, "npc") - unit(0.6, "mm"),
    just = c("left", "top"), gp = title_gp
  )
  popViewport()

  print(p2c, newpage = FALSE, vp = viewport(layout.pos.row = 8, layout.pos.col = 1))
  popViewport()
}

output_suffix <- Sys.getenv("FIG2_OUTPUT_SUFFIX", "_subfigure_titles_bold_ac_trimmed_blue_ac")
output_stem <- paste0("Fig_2", output_suffix)
base <- file.path(out_root, output_stem)
save_pub_grid(draw_fig2, base, 183, 168)

write.csv(fig2_sum, file.path(out_root, paste0(output_stem, "_summary_plot_data_elbow.csv")), row.names = FALSE)
write.csv(patient_occ_summary, file.path(out_root, paste0(output_stem, "_patient_level_summary_plot_data_elbow.csv")), row.names = FALSE)
write.csv(fig2_con, file.path(out_root, paste0(output_stem, "_contrast_plot_data_elbow.csv")), row.names = FALSE)
write.csv(patient_occ, file.path(out_root, paste0(output_stem, "_patient_level_occupancy_elbow.csv")), row.names = FALSE)
write.csv(fig2_patient_con, file.path(out_root, paste0(output_stem, "_patient_level_contrasts_elbow.csv")), row.names = FALSE)
write.csv(patient_ictal_diff, file.path(out_root, paste0(output_stem, "_patient_level_ictal_shift_elbow.csv")), row.names = FALSE)
write.csv(patient_ictal_summary, file.path(out_root, paste0(output_stem, "_patient_level_ictal_shift_summary_elbow.csv")), row.names = FALSE)
write.csv(case_patient_agreement, file.path(out_root, paste0(output_stem, "_case_patient_agreement_elbow.csv")), row.names = FALSE)

message("Rendered ", output_stem, " to ", out_root)
