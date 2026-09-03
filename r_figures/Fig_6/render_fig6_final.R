#!/usr/bin/env Rscript

# Privacy-preserving public renderer for the current four-panel Fig. 6.
# It uses aggregate estimates in Source_Data.xlsx and never reads or writes
# participant-, seizure- or recording-level data.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/")
workbook_path <- if (length(args) >= 1) args[1] else file.path(repo_root, "Source_Data.xlsx")
out_dir <- if (length(args) >= 2) args[2] else file.path(repo_root, "public_outputs", "Fig_6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_block <- function(range) readxl::read_excel(workbook_path, sheet = "Fig6", range = range, col_names = TRUE)
a <- read_block("A3:K27")
b <- read_block("A31:K63")
c <- read_block("A67:K79")
d <- read_block("A83:K91")

stage_levels <- c("Early", "Mid", "Late", "Post")
stage_cols <- c(Early = "#6EAD67", Mid = "#E6A34A", Late = "#D2635E", Post = "#8D7BB8")
prep <- function(x) {
  x$stage <- factor(x$stage, levels = stage_levels)
  for (name in intersect(c("value", "ci_low", "ci_high", "raw_p", "adjusted_p", "n_patients", "B"), names(x))) {
    x[[name]] <- suppressWarnings(as.numeric(x[[name]]))
  }
  x
}
a <- prep(a); b <- prep(b); c <- prep(c); d <- prep(d)
a$metric <- factor(a$metric, levels = rev(c("S", "Gini gap", "Entropy gap", "PR gap", "Reciprocal family", "Dense subset")))

theme_public <- function() {
  theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      strip.background = element_blank(), strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "#E8ECEF"),
      legend.position = "top", legend.title = element_blank()
    )
}

p_a <- ggplot(a, aes(value, metric, colour = stage)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#858B90") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.55) +
  geom_point(aes(fill = adjusted_p < 0.05), shape = 21, size = 2.2, stroke = 0.65) +
  facet_grid(. ~ stage) + scale_colour_manual(values = stage_cols, guide = "none") +
  scale_fill_manual(values = c(`TRUE` = "#2B7A78", `FALSE` = "white"), guide = "none") +
  labs(title = "(a) Density and order matching", x = "Matched difference", y = NULL) + theme_public()

p_b <- ggplot(b, aes(value, metric, colour = stage)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#858B90") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.55) +
  geom_point(shape = 21, fill = "white", size = 2.2, stroke = 0.65) +
  facet_grid(. ~ group) + scale_colour_manual(values = stage_cols) +
  labs(title = "(b) DMWA randomization controls", x = "Empirical minus null", y = NULL) + theme_public()

p_c <- ggplot(c, aes(stage, value, colour = group, group = group)) +
  geom_hline(yintercept = 0, linetype = "22", colour = "#858B90") + geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.08, linewidth = 0.55) +
  geom_point(shape = 21, fill = "white", size = 2.3, stroke = 0.7) +
  scale_colour_manual(values = c(Original = "#2B7A78", `Coverage only` = "#8D7BB8", `degree-bias only` = "#E6A34A")) +
  labs(title = "(c) S component ablation", x = NULL, y = "Change in S") + theme_public()

p_d <- ggplot(d, aes(value, stage, colour = group)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "#858B90") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, linewidth = 0.55) +
  geom_point(shape = 23, fill = "white", size = 2.4, stroke = 0.7) +
  scale_colour_manual(values = c(`Coverage only` = "#8D7BB8", `degree-bias only` = "#E6A34A")) +
  labs(title = "(d) Residual effects", x = "Residual effect", y = NULL) + theme_public()

figure <- p_a / p_b / (p_c | p_d) + plot_layout(heights = c(1.05, 1.0, 1.0))
ggsave(file.path(out_dir, "Fig_6_public_aggregate.pdf"), figure, width = 183, height = 188, units = "mm", device = cairo_pdf)
ggsave(file.path(out_dir, "Fig_6_public_aggregate.png"), figure, width = 183, height = 188, units = "mm", dpi = 300, bg = "white")
message("Created aggregate-only public Fig. 6 outputs in: ", out_dir)
