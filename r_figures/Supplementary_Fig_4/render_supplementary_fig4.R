#!/usr/bin/env Rscript

# Aggregate-only public renderer for Supplementary Fig. 4a-b.

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
out_dir <- if (length(args) >= 2) args[2] else file.path(repo_root, "public_outputs", "Supplementary_Fig_4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- readxl::read_excel(workbook_path, sheet = "SuppFig4", range = "A3:I13", col_names = TRUE)
for (name in c("mean_excess_pp", "ci_low", "ci_high", "adjusted_p")) d[[name]] <- as.numeric(d[[name]])
d$clinical_class <- factor(d$clinical_class, levels = rev(c("NIZ", "PZ", "PZ+NIZ", "mixed SOZ-inclusive", "SOZ-centered")))

make_panel <- function(role_name, colour, label) {
  x <- d[d$role == role_name, ]
  ggplot(x, aes(mean_excess_pp, clinical_class)) +
    geom_vline(xintercept = 0, linetype = "22", colour = "#858B90") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, colour = colour, linewidth = 0.6) +
    geom_point(shape = 21, fill = "white", colour = colour, size = 2.5, stroke = 0.75) +
    geom_text(data = x[x$adjusted_p < 0.05, ], aes(x = ci_high, label = "*"), hjust = -0.45, colour = colour, size = 4) +
    labs(title = label, x = "Role-load share excess (percentage points)", y = NULL) +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0), panel.grid.major.y = element_line(colour = "#E8ECEF", linewidth = 0.2))
}

figure <- make_panel("Source", "#3C6E91", "(a) Source-like load enrichment") |
  make_panel("Sink", "#8A5A7B", "(b) Sink-like load enrichment")
ggsave(file.path(out_dir, "Supplementary_Fig_4_public_aggregate.pdf"), figure, width = 183, height = 82, units = "mm", device = cairo_pdf)
ggsave(file.path(out_dir, "Supplementary_Fig_4_public_aggregate.png"), figure, width = 183, height = 82, units = "mm", dpi = 300, bg = "white")
message("Created aggregate-only Supplementary Fig. 4 outputs in: ", out_dir)
