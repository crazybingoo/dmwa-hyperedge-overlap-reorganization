#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
project_root <- Sys.getenv("PAPER1_PROJECT_ROOT", unset = repo_root)
stats_dir <- file.path(project_root, "stats", "supplementary_sensitivity")
source(file.path(repo_root, "r_figures", "shared", "paper1_figure_style_config.R"))
theme_set(theme_paper1())

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required table: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

save_pub <- function(plot, out_dir, stem, width_mm, height_mm, dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, stem); w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white"); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = dpi, background = "white"); print(plot); dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = dpi, background = "white", compression = "lzw"); print(plot); dev.off()
  data.frame(file = c(paste0(base, c(".svg", ".pdf", ".png", ".tiff"))),
             width_mm = width_mm, height_mm = height_mm, dpi = c(NA, NA, dpi, dpi))
}

exact_signflip_vector <- function(d) {
  n <- nrow(d); observed <- sqrt(sum(colMeans(d)^2)); total <- 2^n
  stats <- numeric(total)
  for (k in 0:(total - 1)) {
    signs <- ifelse(bitwAnd(k, bitwShiftL(1L, 0:(n - 1))) > 0, 1, -1)
    stats[k + 1] <- sqrt(sum(colMeans(d * signs)^2))
  }
  sum(stats >= observed - 1e-12) / total
}

exact_signflip_scalar <- function(x) {
  n <- length(x); observed <- abs(mean(x)); total <- 2^n
  stats <- numeric(total)
  for (k in 0:(total - 1)) {
    signs <- ifelse(bitwAnd(k, bitwShiftL(1L, 0:(n - 1))) > 0, 1, -1)
    stats[k + 1] <- abs(mean(x * signs))
  }
  (sum(stats >= observed - 1e-12) + 1) / (total + 1)
}

boot_vector_norm <- function(d, B = 5000, seed = 20260714) {
  set.seed(seed); n <- nrow(d)
  vals <- replicate(B, {z <- d[sample.int(n, n, replace = TRUE), , drop = FALSE]; sqrt(sum(colMeans(z)^2))})
  unname(quantile(vals, c(0.025, 0.975), na.rm = TRUE))
}

boot_mean <- function(x, B = 5000, seed = 20260714) {
  set.seed(seed); n <- length(x)
  vals <- replicate(B, mean(x[sample.int(n, n, replace = TRUE)], na.rm = TRUE))
  unname(quantile(vals, c(0.025, 0.975), na.rm = TRUE))
}

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

# Supplementary Fig. 1 -----------------------------------------------------
win1 <- read_required(file.path(stats_dir, "supp_fig1_synchrony_sensitivity_window.csv"))
rules <- read_required(file.path(stats_dir, "supp_fig1_case_rules.csv"))
features <- c("S", "source_sink_gini_gap", "sink_minus_source_entropy_gap", "sink_minus_source_pr_gap")

patient_stage <- win1 %>%
  group_by(variant, case_id, patient_id, stage) %>%
  summarise(across(all_of(features), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  group_by(variant, patient_id, stage) %>%
  summarise(across(all_of(features), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  group_by(variant) %>%
  mutate(across(all_of(features), ~as.numeric(scale(.x)))) %>%
  ungroup()

pre <- patient_stage %>% filter(stage == "pre-ictal") %>% select(variant, patient_id, all_of(features))
ict <- patient_stage %>% filter(stage %in% c("early", "mid", "late")) %>%
  group_by(variant, patient_id) %>% summarise(across(all_of(features), ~mean(.x, na.rm = TRUE)), .groups = "drop")
paired <- inner_join(ict, pre, by = c("variant", "patient_id"), suffix = c("_ictal", "_pre"))
for (f in features) paired[[paste0("d_", f)]] <- paired[[paste0(f, "_ictal")]] - paired[[paste0(f, "_pre")]]

variant_order <- c("PLV baseline", "PLV pre-ictal threshold", "PLV threshold -10%", "PLV threshold +10%",
                   "PLV fixed density", "PLI fixed density", "wPLI fixed density")
variant_labels <- c("Baseline PLV", "Pre-ictal PLV threshold", "PLV threshold -10%", "PLV threshold +10%",
                    "PLV at fixed density", "PLI at fixed density", "wPLI at fixed density")

fig1_summary <- lapply(seq_along(variant_order), function(i) {
  z <- paired %>% filter(variant == variant_order[i])
  d <- as.matrix(z[, paste0("d_", features)])
  ci <- boot_vector_norm(d, seed = 20260714 + i)
  data.frame(variant = variant_order[i], label = variant_labels[i], n_patients = nrow(d),
             centroid_distance = sqrt(sum(colMeans(d)^2)), ci_low = ci[1], ci_high = ci[2],
             exact_p = exact_signflip_vector(d))
}) %>%
  bind_rows() %>%
  mutate(holm_adjusted_p = p.adjust(exact_p, method = "holm"))

component_summary <- paired %>%
  select(variant, patient_id, starts_with("d_")) %>%
  pivot_longer(starts_with("d_"), names_to = "feature", values_to = "difference") %>%
  group_by(variant, feature) %>%
  summarise(
    n_patients = n(),
    mean_difference = mean(difference),
    values = list(difference),
    .groups = "drop"
  ) %>%
  mutate(
    seed = 20260714 + 10L * match(variant, variant_order) + match(feature, paste0("d_", features))
  ) %>%
  rowwise() %>%
  mutate(
    bootstrap_ci = list(boot_mean(unlist(values), seed = seed)),
    ci_low = bootstrap_ci[[1]],
    ci_high = bootstrap_ci[[2]]
  ) %>%
  ungroup() %>%
  select(-values, -seed, -bootstrap_ci) %>%
  mutate(feature = recode(feature,
    d_S = "S",
    d_source_sink_gini_gap = "Source-sink Gini gap",
    d_sink_minus_source_entropy_gap = "Sink-source entropy gap",
    d_sink_minus_source_pr_gap = "Sink-source PR gap"))

write.csv(fig1_summary, file.path(stats_dir, "supp_fig1_synchrony_sensitivity_patient_summary.csv"), row.names = FALSE)
write.csv(component_summary, file.path(stats_dir, "supp_fig1_synchrony_sensitivity_component_shifts.csv"), row.names = FALSE)

rules_long <- rules %>%
  select(case_id, baseline_threshold, preictal_only_threshold) %>%
  pivot_longer(-case_id, names_to = "rule", values_to = "threshold") %>%
  mutate(rule = recode(rule, baseline_threshold = "All-stage", preictal_only_threshold = "Pre-ictal only"),
         rule = factor(rule, levels = c("All-stage", "Pre-ictal only")))

p1a <- ggplot(rules_long, aes(rule, threshold, group = case_id)) +
  geom_line(colour = paper1_alpha(sig("neutral_mid"), 0.42), linewidth = 0.35) +
  geom_point(shape = 21, size = 2.0, stroke = 0.55, fill = "white", colour = sig("neutral_dark")) +
  stat_summary(aes(group = 1), fun = median, geom = "crossbar", width = 0.48,
               colour = sig("teal"), linewidth = 0.55) +
  labs(title = NULL, x = NULL, y = expression(theta[q]^"*")) +
  coord_cartesian(ylim = range(rules_long$threshold) + c(-0.02, 0.02))

p1b <- fig1_summary %>% mutate(label = factor(label, levels = rev(variant_labels))) %>%
  ggplot(aes(centroid_distance, label)) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.55, colour = sig("neutral_mid")) +
  geom_point(shape = 21, size = 2.8, stroke = 0.65, fill = "white", colour = sig("teal")) +
  labs(title = NULL, x = expression(D[s]~"(standardized units)"), y = NULL) +
  theme(panel.grid.major.y = element_blank())

p1c <- component_summary %>%
  mutate(
    variant = factor(variant, levels = variant_order),
    feature = factor(
      feature,
      levels = c("S", "Source-sink Gini gap", "Sink-source entropy gap", "Sink-source PR gap"),
      labels = c("S", "D[G]", "D[H]", "D[PR]")
    )
  ) %>%
  ggplot(aes(mean_difference, variant)) +
  geom_segment(
    aes(x = 0, xend = mean_difference, yend = variant),
    linewidth = 0.34,
    colour = paper1_alpha(sig("neutral_mid"), 0.58)
  ) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0,
    linewidth = 0.38,
    colour = sig("neutral_mid")
  ) +
  geom_point(
    shape = 21, size = 2.15, stroke = 0.58,
    fill = "white", colour = sig("teal")
  ) +
  geom_text(
    aes(label = sprintf("%.2f", mean_difference)),
    nudge_y = 0.32,
    hjust = 0.5,
    size = FIG_HEATMAP_TEXT_SIZE,
    colour = sig("neutral_dark")
  ) +
  facet_grid(. ~ feature, labeller = label_parsed) +
  scale_x_continuous(limits = c(-0.20, 2.30), breaks = c(0, 0.5, 1.0, 1.5, 2.0), expand = expansion(mult = c(0, 0.01))) +
  scale_y_discrete(labels = setNames(variant_labels, variant_order)) +
  labs(title = NULL, x = expression(bar(Delta*z)[list(s,k)]), y = NULL) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = FIG_TEXT_PT),
    panel.grid.major.x = element_line(colour = paper1_alpha(sig("neutral_light"), 0.55), linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.x = grid::unit(8, "pt")
  )

title_strip <- function(label) {
  wrap_elements(
    full = grid::textGrob(
      label,
      x = grid::unit(3.0, "mm"),
      y = grid::unit(0.50, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(
        fontfamily = "Arial",
        fontsize = FIG_PANEL_PT,
        fontface = "bold",
        col = sig("neutral_dark")
      )
    ),
    clip = FALSE
  )
}

panel_a <- title_strip("(a) Threshold rule comparison") / p1a + plot_layout(heights = c(0.13, 1))
panel_b <- title_strip("(b) Reconstruction sensitivity") / p1b + plot_layout(heights = c(0.13, 1))
panel_c <- title_strip("(c) Direction of component shifts") / p1c + plot_layout(heights = c(0.11, 1))

top1 <- wrap_plots(
  wrap_elements(full = panel_a, clip = FALSE),
  wrap_elements(full = panel_b, clip = FALSE),
  ncol = 2,
  widths = c(0.94, 1.06)
)
fig1 <- wrap_plots(
  wrap_elements(full = top1, clip = FALSE),
  wrap_elements(full = panel_c, clip = FALSE),
  ncol = 1,
  heights = c(1, 1.18)
)
out1 <- file.path(project_root, "nature_fig", "Supplementary_Fig_1")
manifest1 <- save_pub(fig1, out1, "Supplementary_Fig_1", 183, 150)
write.csv(manifest1, file.path(out1, "Supplementary_Fig_1_export_manifest.csv"), row.names = FALSE)

if ("--fig1-only" %in% commandArgs(trailingOnly = TRUE)) {
  cat("Supplementary Fig. 1 rendered.\n")
  print(fig1_summary)
  quit(save = "no", status = 0)
}

# Supplementary Fig. 2 -----------------------------------------------------
win2 <- read_required(file.path(stats_dir, "supp_fig2_reciprocity_ratio_window.csv"))
motif_features <- c("simple_fraction", "reciprocal_fraction", "dense_fraction")
patient_stage2 <- win2 %>%
  group_by(reciprocity_ratio, case_id, patient_id, stage) %>%
  summarise(across(all_of(motif_features), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  group_by(reciprocity_ratio, patient_id, stage) %>%
  summarise(across(all_of(motif_features), ~mean(.x, na.rm = TRUE)), .groups = "drop")
pre2 <- patient_stage2 %>% filter(stage == "pre-ictal") %>% select(reciprocity_ratio, patient_id, all_of(motif_features))
ict2 <- patient_stage2 %>% filter(stage %in% c("early", "mid", "late")) %>%
  group_by(reciprocity_ratio, patient_id) %>% summarise(across(all_of(motif_features), ~mean(.x, na.rm = TRUE)), .groups = "drop")
paired2 <- inner_join(ict2, pre2, by = c("reciprocity_ratio", "patient_id"), suffix = c("_ictal", "_pre"))
for (f in motif_features) paired2[[paste0("d_", f)]] <- paired2[[paste0(f, "_ictal")]] - paired2[[paste0(f, "_pre")]]

long2 <- paired2 %>% select(reciprocity_ratio, patient_id, starts_with("d_")) %>%
  pivot_longer(starts_with("d_"), names_to = "metric", values_to = "difference") %>%
  mutate(metric = recode(metric, d_simple_fraction = "Simple", d_reciprocal_fraction = "Reciprocal", d_dense_fraction = "Dense"))
sum2 <- long2 %>% group_by(reciprocity_ratio, metric) %>% group_modify(~{
  ci <- boot_mean(.x$difference, seed = 20260714 + round(.y$reciprocity_ratio * 100) + match(.y$metric, c("Simple","Reciprocal","Dense")))
  data.frame(n_patients = length(.x$difference), mean_difference = mean(.x$difference), ci_low = ci[1], ci_high = ci[2],
             exact_p = exact_signflip_scalar(.x$difference))
}) %>% ungroup() %>% group_by(metric) %>% mutate(adjusted_p = p.adjust(exact_p, method = "BH")) %>% ungroup()
write.csv(sum2, file.path(stats_dir, "supp_fig2_reciprocity_ratio_patient_summary.csv"), row.names = FALSE)

make_ratio_panel <- function(metric_name, panel, colour) {
  raw <- long2 %>% filter(metric == metric_name)
  sm <- sum2 %>% filter(metric == metric_name)
  ggplot(raw, aes(reciprocity_ratio, difference, group = patient_id)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = sig("neutral_mid")) +
    geom_vline(xintercept = 0.85, linewidth = 0.35, linetype = 2, colour = sig("neutral_mid")) +
    geom_line(colour = paper1_alpha(colour, 0.24), linewidth = 0.32) +
    geom_point(shape = 21, size = 1.15, stroke = 0.3, fill = "white", colour = paper1_alpha(colour, 0.48)) +
    geom_errorbar(data = sm, aes(y = mean_difference, ymin = ci_low, ymax = ci_high, group = NULL),
                  width = 0.012, linewidth = 0.55, colour = sig("neutral_dark")) +
    geom_line(data = sm, aes(y = mean_difference, group = 1), linewidth = 0.7, colour = colour) +
    geom_point(data = sm, aes(y = mean_difference, group = NULL), shape = 21, size = 2.2, stroke = 0.55,
               fill = "white", colour = colour) +
    scale_x_continuous(breaks = seq(0.70, 0.95, 0.05)) +
    labs(title = NULL, x = expression(tau[rho]), y = expression(Delta*f[F])) +
    theme(panel.grid.major.x = element_blank())
}

p2a <- make_ratio_panel("Simple", "a", sig("blue"))
p2b <- make_ratio_panel("Reciprocal", "b", sig("red"))
p2c <- make_ratio_panel("Dense", "c", sig("purple"))
panel2a <- title_strip("(a) Simple motifs") / p2a + plot_layout(heights = c(0.13, 1))
panel2b <- title_strip("(b) Reciprocal motifs") / p2b + plot_layout(heights = c(0.13, 1))
panel2c <- title_strip("(c) Dense motifs") / p2c + plot_layout(heights = c(0.13, 1))
fig2 <- wrap_plots(
  wrap_elements(full = panel2a, clip = FALSE),
  wrap_elements(full = panel2b, clip = FALSE),
  wrap_elements(full = panel2c, clip = FALSE),
  ncol = 3
)
out2 <- file.path(project_root, "nature_fig", "Supplementary_Fig_2")
manifest2 <- save_pub(fig2, out2, "Supplementary_Fig_2", 183, 82)
write.csv(manifest2, file.path(out2, "Supplementary_Fig_2_export_manifest.csv"), row.names = FALSE)

cat("Supplementary Figs. 1 and 2 rendered.\n")
print(fig1_summary)
print(sum2)
