# Shared visual system for paper 1 Fig. 2-8.
# This file only controls colours, typography, line weights and export-facing theme defaults.
# It does not change data processing, statistics, coordinates or scientific conclusions.

FIG_TEXT_PT <- 9.8
FIG_PANEL_PT <- 11.8
FIG_HEATMAP_TEXT_PT <- 7.6
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276
FIG_HEATMAP_TEXT_SIZE <- FIG_HEATMAP_TEXT_PT / 2.845276

paper1_palette <- list(
  primary_teal = "#2B7A78",
  light_teal_fill = "#2B7A78",
  accent_red = "#D2635E",
  background = "#FFFFFF",
  text_dark = "#222A33",
  axis_dark = "#30343B",
  grid_light = "#EDF3F3",
  auxiliary_gray = "#CFD2D2",
  annotation_gray = "#5A5F66",
  pre = "#5D83B5",
  early = "#6EAD67",
  mid = "#E6A34A",
  late = "#D2635E",
  post = "#8D7BB8"
)

paper1_alpha <- function(x, alpha = 0.2) {
  grDevices::adjustcolor(x, alpha.f = alpha)
}

stage_palette <- c(
  "pre-ictal" = paper1_palette$pre,
  "early" = paper1_palette$early,
  "mid" = paper1_palette$mid,
  "late" = paper1_palette$late,
  "post-ictal" = paper1_palette$post,
  "ictal-combined" = paper1_palette$primary_teal
)

signal_palette <- c(
  neutral_dark = paper1_palette$text_dark,
  neutral_mid = paper1_palette$annotation_gray,
  neutral_light = paper1_palette$auxiliary_gray,
  charcoal = paper1_palette$text_dark,
  grey = paper1_palette$annotation_gray,
  light_grey = paper1_palette$auxiliary_gray,
  blue = paper1_palette$pre,
  teal = paper1_palette$primary_teal,
  red = paper1_palette$accent_red,
  orange = paper1_palette$mid,
  green = paper1_palette$early,
  purple = paper1_palette$post,
  simple = paper1_palette$pre,
  ffl = paper1_palette$mid,
  reciprocal = paper1_palette$accent_red,
  dense = paper1_palette$post,
  neutral = paper1_palette$text_dark,
  primary_teal = paper1_palette$primary_teal,
  accent_red = paper1_palette$accent_red,
  auxiliary_gray = paper1_palette$auxiliary_gray
)

sig <- function(name) unname(signal_palette[name])

feature_set_palette <- c(
  "undirected_plus_order" = paper1_palette$auxiliary_gray,
  "undirected_linegraph" = paper1_palette$post,
  "hyperedge_density_order" = paper1_palette$mid,
  "DMWA_direction_role" = paper1_palette$primary_teal,
  "pairwise_plv_density" = paper1_palette$pre
)

model_palette <- c(
  "role_reciprocity_hybrid" = paper1_palette$mid,
  "density_only" = paper1_palette$pre,
  "source_concentration" = paper1_palette$primary_teal,
  "overlap_synchronization" = paper1_palette$early,
  "weight_strengthening" = paper1_palette$annotation_gray,
  "reciprocity_only" = paper1_palette$post
)

stage_fill_palette <- paper1_alpha(stage_palette, 0.18)
feature_set_fill_palette <- paper1_alpha(feature_set_palette, 0.22)
model_fill_palette <- paper1_alpha(model_palette, 0.14)

diverging_low <- "#8FA8C8"
diverging_mid <- "#FAFAF8"
diverging_high <- paper1_palette$accent_red

theme_paper1 <- function(base_size = FIG_TEXT_PT, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.background = element_rect(fill = paper1_palette$background, colour = NA),
      panel.background = element_rect(fill = paper1_palette$background, colour = NA),
      axis.line = element_line(linewidth = 0.36, colour = paper1_palette$axis_dark),
      axis.ticks = element_line(linewidth = 0.32, colour = paper1_palette$axis_dark),
      axis.title = element_text(size = base_size, colour = paper1_palette$text_dark),
      axis.title.x = element_text(margin = margin(t = 2.0)),
      axis.title.y = element_text(margin = margin(r = 1.5)),
      axis.text = element_text(size = base_size, colour = paper1_palette$text_dark),
      legend.title = element_text(size = base_size, colour = paper1_palette$text_dark),
      legend.text = element_text(size = base_size, colour = paper1_palette$text_dark),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = paper1_palette$text_dark),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0, colour = paper1_palette$text_dark),
      plot.subtitle = element_text(size = base_size, colour = paper1_palette$annotation_gray),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = paper1_palette$grid_light),
      panel.grid.major.x = element_line(linewidth = 0.16, colour = paper1_palette$grid_light),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

theme_set(theme_paper1())
panel_tag_theme <- theme(plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", colour = paper1_palette$text_dark))
