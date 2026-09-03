library(ggplot2)
library(dplyr)
library(patchwork)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
project_root <- Sys.getenv("PAPER1_PROJECT_ROOT", unset = repo_root)
stats_dir <- file.path(project_root, "stats")
out_dir <- file.path(project_root, "nature_fig", "Fig_5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "r_figures", "shared", "paper1_figure_style_config.R"))

read_stats <- function(filename) {
  read.csv(
    file.path(stats_dir, filename),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

stage_levels <- c("Early", "Mid", "Late", "Post", "Combined")
stage_map <- c(
  "early" = "Early",
  "mid" = "Mid",
  "late" = "Late",
  "post-ictal" = "Post",
  "ictal-combined" = "Combined"
)
stage_map_with_pre <- c("pre-ictal" = "Pre", stage_map)

feature_order <- c("S", "D_G", "D_H", "D_PR", "f_recip", "f_dense")
feature_map <- c(
  "S" = "S",
  "source_sink_gini_gap" = "D_G",
  "sink_minus_source_entropy_gap" = "D_H",
  "sink_minus_source_pr_gap" = "D_PR",
  "reciprocal_fraction_mean" = "f_recip",
  "dense_fraction_mean" = "f_dense"
)
feature_family_map <- c(
  "S" = "Symmetry",
  "D_G" = "Role",
  "D_H" = "Role",
  "D_PR" = "Role",
  "f_recip" = "Motif",
  "f_dense" = "Motif"
)

centroid <- read_stats(
  "phase17_elbow_joint_state_space_centroid_tests.csv"
) %>%
  mutate(
    stage_label = factor(unname(stage_map[target_stage]), levels = stage_levels),
    stage_y = as.numeric(factor(stage_label, levels = rev(stage_levels))),
    value_label = sprintf("%.2f", centroid_distance_z),
    label_x = bootstrap_ci_high_z + 0.08
  )

feature_shift <- read_stats(
  "phase17_elbow_joint_state_space_feature_shifts.csv"
) %>%
  filter(target_stage == "ictal-combined") %>%
  mutate(
    feature_label = factor(
      unname(feature_map[feature]),
      levels = rev(feature_order)
    ),
    feature_family = factor(
      unname(feature_family_map[as.character(feature_label)]),
      levels = c("Symmetry", "Role", "Motif")
    ),
    value_label = sprintf("%+.2f", mean_diff_z),
    label_x = mean_diff_z + 0.07,
    label_hjust = 0
  )

feature_set_map <- c(
  "DMWA_full" = "Six-feature DMWA",
  "block_balanced" = "Block-balanced",
  "S_only" = "S only",
  "role_only" = "Role only",
  "motif_only" = "Motif only",
  "density_order_only" = "Density/order"
)
feature_set_levels <- c(
  "Six-feature DMWA",
  "Block-balanced",
  "S only",
  "Role only",
  "Motif only",
  "Density/order"
)

feature_sets <- read_stats("results_logic_feature_set_rms.csv") %>%
  filter(target_stage == "ictal-combined") %>%
  mutate(
    feature_set_label = factor(
      unname(feature_set_map[feature_set]),
      levels = rev(feature_set_levels)
    ),
    feature_set_y = as.numeric(feature_set_label),
    rms_ci_low_z = bootstrap_ci_low_z / sqrt(n_features),
    rms_ci_high_z = bootstrap_ci_high_z / sqrt(n_features),
    value_label = sprintf("%.2f", rms_distance_z),
    label_x = rms_ci_high_z + 0.04
  )

state_features <- names(feature_map)
state_z_columns <- paste0(state_features, "_z")
state_points <- read_stats(
  "phase17_elbow_joint_state_space_patient_stage_features.csv"
) %>%
  filter(stage %in% names(stage_map_with_pre)) %>%
  select(patient_id, stage, all_of(state_z_columns))
names(state_points)[match(state_z_columns, names(state_points))] <- state_features
state_points <- state_points %>%
  mutate(
    stage_label = factor(
      unname(stage_map_with_pre[stage]),
      levels = c("Pre", "Early", "Mid", "Late", "Post")
    )
  )

pre_state <- state_points %>%
  filter(stage_label == "Pre") %>%
  select(patient_id, all_of(state_features))

state_distance <- state_points %>%
  filter(stage_label != "Pre") %>%
  left_join(pre_state, by = "patient_id", suffix = c("", "_pre")) %>%
  rowwise() %>%
  mutate(
    individual_distance_z = sqrt(sum(
      (
        c_across(all_of(state_features)) -
          c_across(all_of(paste0(state_features, "_pre")))
      )^2
    ))
  ) %>%
  ungroup() %>%
  mutate(
    stage_label = factor(
      stage_label,
      levels = c("Early", "Mid", "Late", "Post")
    )
  )

state_distance_summary <- state_distance %>%
  group_by(stage_label) %>%
  summarise(
    n_patients = n_distinct(patient_id),
    min_distance_z = min(individual_distance_z),
    q1_distance_z = quantile(individual_distance_z, 0.25),
    median_distance_z = median(individual_distance_z),
    q3_distance_z = quantile(individual_distance_z, 0.75),
    max_distance_z = max(individual_distance_z),
    .groups = "drop"
  )

rain_levels <- c("Post", "Late", "Mid", "Early")
state_distance_plot <- state_distance %>%
  mutate(
    stage_label_plot = factor(as.character(stage_label), levels = rain_levels),
    stage_y = as.numeric(stage_label_plot)
  )
state_distance_summary_plot <- state_distance_summary %>%
  mutate(
    stage_label_plot = factor(as.character(stage_label), levels = rain_levels),
    stage_y = as.numeric(stage_label_plot)
  )

make_half_density <- function(
  data,
  value_col = "individual_distance_z",
  group_col = "stage_label_plot",
  y_col = "stage_y",
  width = 0.30
) {
  x_max <- max(data[[value_col]], na.rm = TRUE) * 1.04
  pieces <- lapply(levels(data[[group_col]]), function(stage_name) {
    stage_data <- data[data[[group_col]] == stage_name, , drop = FALSE]
    values <- stage_data[[value_col]]
    values <- values[is.finite(values)]
    density <- stats::density(
      values,
      from = 0,
      to = x_max,
      n = 180,
      adjust = 0.95
    )
    y_base <- unique(stage_data[[y_col]])[1]
    y_top <- y_base + density$y / max(density$y) * width
    rbind(
      data.frame(stage_label_plot = stage_name, x = density$x, y = y_top),
      data.frame(
        stage_label_plot = stage_name,
        x = rev(density$x),
        y = rep(y_base, length(density$x))
      )
    )
  })
  bind_rows(pieces) %>%
    mutate(
      stage_label_plot = factor(stage_label_plot, levels = levels(data[[group_col]]))
    )
}

state_raincloud <- make_half_density(state_distance_plot)

stage_colours <- c(
  "Pre" = "#8C8C8C",
  "Early" = "#4C72B0",
  "Mid" = "#5A9E90",
  "Late" = "#C44E52",
  "Post" = "#8172B2",
  "Combined" = "#262626"
)
family_colours <- c(
  "Symmetry" = "#4C72B0",
  "Role" = "#5A9E90",
  "Motif" = "#C44E52"
)
set_colours <- c(
  "Six-feature DMWA" = "#262626",
  "Block-balanced" = "#7A5C8E",
  "S only" = "#4C72B0",
  "Role only" = "#5A9E90",
  "Motif only" = "#C44E52",
  "Density/order" = "#B9A549"
)

panel_title_theme <- theme(
  plot.title.position = "plot",
  plot.title = element_text(
    size = FIG_PANEL_PT,
    face = "bold",
    hjust = 0,
    colour = paper1_palette$text_dark
  )
)
lucas_theme <- theme(
  panel.grid = element_blank(),
  axis.line = element_line(linewidth = 0.36, colour = "#30343B"),
  axis.ticks = element_line(linewidth = 0.32, colour = "#30343B"),
  legend.position = "none"
)

p_a <- ggplot() +
  geom_polygon(
    data = state_raincloud,
    aes(x = x, y = y, group = stage_label_plot, fill = stage_label_plot),
    alpha = 0.16,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_segment(
    data = state_distance_summary_plot,
    aes(
      x = min_distance_z,
      xend = max_distance_z,
      y = stage_y,
      yend = stage_y
    ),
    colour = "#404040",
    linewidth = 0.34
  ) +
  geom_rect(
    data = state_distance_summary_plot,
    aes(
      xmin = q1_distance_z,
      xmax = q3_distance_z,
      ymin = stage_y - 0.10,
      ymax = stage_y + 0.10
    ),
    fill = "white",
    colour = "#404040",
    linewidth = 0.40
  ) +
  geom_segment(
    data = state_distance_summary_plot,
    aes(
      x = median_distance_z,
      xend = median_distance_z,
      y = stage_y - 0.12,
      yend = stage_y + 0.12
    ),
    colour = "#404040",
    linewidth = 0.58
  ) +
  geom_point(
    data = state_distance_plot,
    aes(
      x = individual_distance_z,
      y = stage_y - 0.24,
      colour = stage_label_plot
    ),
    position = position_jitter(width = 0, height = 0.055, seed = 20260717),
    size = 1.48,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = paper1_alpha(stage_colours[rain_levels], 0.24)) +
  scale_colour_manual(values = stage_colours) +
  scale_x_continuous(
    limits = c(0, max(state_distance$individual_distance_z) * 1.10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq_along(rain_levels),
    labels = rain_levels,
    limits = c(0.55, length(rain_levels) + 0.55)
  ) +
  labs(
    title = "(a) Patient-level displacement",
    x = expression(delta[list(p,s)]~"(z)"),
    y = NULL
  ) +
  theme_paper1() +
  panel_title_theme +
  lucas_theme +
  theme(plot.margin = margin(4, 10, 6, 4))

component_lower <- min(0, feature_shift$mean_diff_z) - 0.15
p_b <- ggplot(feature_shift, aes(mean_diff_z, feature_label)) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.30,
    linetype = "dashed",
    colour = "#B0B0B0"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = mean_diff_z,
      y = feature_label,
      yend = feature_label,
      colour = feature_family
    ),
    linewidth = 0.58,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = feature_family),
    shape = 21,
    fill = "white",
    size = 2.45,
    stroke = 0.68,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = label_x, label = value_label, hjust = label_hjust),
    size = FIG_HEATMAP_TEXT_SIZE,
    family = "Arial",
    colour = paper1_palette$text_dark
  ) +
  scale_colour_manual(values = family_colours) +
  scale_x_continuous(
    limits = c(component_lower, max(feature_shift$mean_diff_z) * 1.30),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(
    labels = parse(text = c("f[dense]", "f[recip]", "D[PR]", "D[H]", "D[G]", "S"))
  ) +
  labs(
    title = "(b) Combined-ictal component shifts",
    x = expression(bar(Delta*z)[list(s,k)]),
    y = NULL
  ) +
  theme_paper1() +
  panel_title_theme +
  lucas_theme

p_c <- ggplot(centroid, aes(centroid_distance_z, stage_y)) +
  geom_segment(
    aes(
      x = bootstrap_ci_low_z,
      xend = bootstrap_ci_high_z,
      y = stage_y,
      yend = stage_y,
      colour = stage_label
    ),
    linewidth = 0.54,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = stage_label),
    shape = 21,
    fill = "white",
    size = 2.45,
    stroke = 0.68,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = label_x, label = value_label),
    hjust = 0,
    size = FIG_HEATMAP_TEXT_SIZE,
    family = "Arial",
    colour = paper1_palette$text_dark
  ) +
  scale_colour_manual(values = stage_colours) +
  scale_x_continuous(
    limits = c(0, max(centroid$bootstrap_ci_high_z) * 1.27),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq_along(rev(stage_levels)),
    labels = rev(stage_levels),
    limits = c(0.65, length(stage_levels) + 0.45)
  ) +
  labs(
    title = "(c) Group centroid displacement",
    x = expression(D[s]~"(z)"),
    y = NULL
  ) +
  theme_paper1() +
  panel_title_theme +
  lucas_theme

six_feature_rms <- feature_sets$rms_distance_z[
  feature_sets$feature_set_label == "Six-feature DMWA"
]
p_d <- ggplot(
  feature_sets,
  aes(rms_distance_z, feature_set_y, colour = feature_set_label)
) +
  geom_vline(
    xintercept = six_feature_rms,
    linewidth = 0.30,
    linetype = "dotted",
    colour = "#A6ADB5"
  ) +
  geom_segment(
    aes(
      x = rms_ci_low_z,
      xend = rms_ci_high_z,
      y = feature_set_y,
      yend = feature_set_y
    ),
    linewidth = 0.56,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(
    shape = 21,
    fill = "white",
    size = 2.45,
    stroke = 0.68,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = label_x, label = value_label),
    hjust = 0,
    size = FIG_HEATMAP_TEXT_SIZE,
    family = "Arial",
    colour = paper1_palette$text_dark
  ) +
  scale_colour_manual(values = set_colours) +
  scale_x_continuous(
    limits = c(0, max(feature_sets$rms_ci_high_z) * 1.19),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq_along(rev(feature_set_levels)),
    labels = rev(feature_set_levels),
    limits = c(0.65, length(feature_set_levels) + 0.45)
  ) +
  labs(
    title = "(d) Dimension-normalized sensitivity",
    x = expression(D[s]^{RMS}~"(z)"),
    y = NULL
  ) +
  theme_paper1() +
  panel_title_theme +
  lucas_theme

figure_design <- "
AABB
AACC
AADD
"
fig_5 <- p_a + p_b + p_c + p_d +
  plot_layout(
    design = figure_design,
    widths = c(1.10, 1.10, 1.00, 1.00),
    heights = c(1.00, 0.90, 0.94),
    guides = "keep"
  )

width_mm <- 183
height_mm <- 150
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
base <- file.path(out_dir, "Fig_5")

svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(fig_5)
grDevices::dev.off()

grDevices::cairo_pdf(
  paste0(base, ".pdf"),
  width = width_in,
  height = height_in,
  family = "Arial"
)
print(fig_5)
grDevices::dev.off()

ragg::agg_tiff(
  paste0(base, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  compression = "lzw"
)
print(fig_5)
grDevices::dev.off()

ragg::agg_png(
  paste0(base, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300
)
print(fig_5)
grDevices::dev.off()

source_a <- state_distance %>%
  transmute(
    panel = "a",
    patient_id = patient_id,
    stage = as.character(stage_label),
    individual_distance_z = individual_distance_z
  )
source_b <- feature_shift %>%
  transmute(
    panel = "b",
    target_stage = target_stage,
    feature = feature,
    mean_diff_z = mean_diff_z,
    mean_diff_raw = mean_diff_raw
  )
source_c <- centroid %>%
  transmute(
    panel = "c",
    target_stage = target_stage,
    n_patients = n_patients,
    centroid_distance_z = centroid_distance_z,
    bootstrap_ci_low_z = bootstrap_ci_low_z,
    bootstrap_ci_high_z = bootstrap_ci_high_z,
    exact_signflip_p = exact_signflip_p,
    secondary_bh_adjusted_p = secondary_bh_adjusted_p
  )
source_d <- feature_sets %>%
  transmute(
    panel = "d",
    feature_set = feature_set,
    n_patients = n_patients,
    n_features = n_features,
    centroid_distance_z = centroid_distance_z,
    rms_distance_z = rms_distance_z,
    rms_ci_low_z = rms_ci_low_z,
    rms_ci_high_z = rms_ci_high_z,
    exact_signflip_p = exact_signflip_p
  )
source_data <- bind_rows(source_a, source_b, source_c, source_d)
write.csv(
  source_data,
  file.path(out_dir, "Fig_5_source_data.csv"),
  row.names = FALSE,
  na = ""
)

statistics_manifest <- data.frame(
  panel = c("a", "b", "c", "d"),
  analysis_unit = c(
    "patient",
    "patient",
    "patient",
    "patient"
  ),
  n = c(14, 14, 14, 14),
  estimand = c(
    "individual Euclidean distance from matched pre-ictal",
    "mean paired standardized coordinate shift, combined ictal vs pre-ictal",
    "Euclidean norm of the group-mean paired difference vector",
    "dimension-normalized RMS of the combined-ictal centroid displacement"
  ),
  uncertainty = c(
    "distribution, median and interquartile range",
    "not displayed",
    "95% patient-bootstrap CI with standardization re-estimated per resample",
    "95% patient-bootstrap CI divided by sqrt(number of coordinates)"
  ),
  inference = c(
    "descriptive",
    "descriptive",
    "exact sign-flip; four stage-specific tests BH-adjusted",
    "sensitivity summaries; not a ranking of representations"
  )
)
write.csv(
  statistics_manifest,
  file.path(out_dir, "Fig_5_statistics_manifest.csv"),
  row.names = FALSE
)

export_manifest <- data.frame(
  file = paste0(base, c(".svg", ".pdf", ".tiff", ".png")),
  width_mm = width_mm,
  height_mm = height_mm,
  dpi = c(NA, NA, 600, 300)
)
write.csv(
  export_manifest,
  file.path(out_dir, "Fig_5_export_manifest.csv"),
  row.names = FALSE
)

contract_lines <- c(
  "# Fig. 5 figure contract",
  "",
  "- Core conclusion: the exploratory six-feature DMWA state shows a coordinated patient-level displacement from pre-ictal baseline.",
  "- Sensitivity boundary: the conclusion persists after the three correlated role metrics are collapsed into one re-standardized role block.",
  "- Comparison boundary: dimension-normalized summaries are sensitivity diagnostics, not evidence that one representation is superior.",
  "- Analysis unit: 14 patients; combined ictal averages early, mid and late within patient.",
  "- Inference: combined ictal is the primary exact sign-flip test; the four stage-specific tests form one Benjamini-Hochberg family.",
  "- Bootstrap: patients are resampled and feature standardization is re-estimated in every resample.",
  "- Typography: FIG_TEXT_PT = 9.8; FIG_PANEL_PT = 11.8; FIG_HEATMAP_TEXT_PT = 7.6.",
  "- Export: 183 x 150 mm; SVG, PDF, 600-dpi TIFF and 300-dpi PNG."
)
writeLines(contract_lines, file.path(out_dir, "Fig_5_figure_contract.md"))

message("Rendered final Fig. 5 and wrote source/QA manifests.")
