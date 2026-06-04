#!/usr/bin/env Rscript
# Render phylogenetic trees for FungiFlow (rectangular or circular).
# Usage: Rscript plot_tree.R <treefile> <output.png> <layout>
#   layout: rectangular | circular  (default: rectangular)

suppressPackageStartupMessages({
  library(ggtree)
  library(ggplot2)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_tree.R <treefile> <output.png> [rectangular|circular]")
}

tree_file <- args[1]
output_image <- args[2]
layout <- if (length(args) >= 3) tolower(args[3]) else "rectangular"

if (!file.exists(tree_file)) {
  stop("Tree file not found: ", tree_file)
}

tree <- read.tree(tree_file)
num_tips <- length(tree$tip.label)

tip_size <- if (num_tips > 80) 1.0 else if (num_tips > 40) 1.3 else 1.8
support_cutoff <- 70

node_support <- function(tree) {
  if (is.null(tree$node.label)) {
    return(rep(NA_real_, Nnode(tree)))
  }
  suppressWarnings(as.numeric(tree$node.label))
}

tree$node.label <- as.character(tree$node.label)

if (layout == "circular") {
  p <- ggtree(tree, layout = "circular", open = TRUE, size = 0.25) +
    geom_tiplab2(
      aes(angle = angle),
      size = tip_size,
      color = "#0a160a",
      align = TRUE,
      linesize = 0.15
    ) +
    geom_nodelab2(
      aes(
        subset = !is.na(label) & label != "" & suppressWarnings(as.numeric(label)) >= support_cutoff,
        label = label
      ),
      size = 1.2,
      color = "#1b4332"
    ) +
    theme(
      plot.margin = margin(8, 8, 8, 8),
      plot.background = element_rect(fill = "white", color = NA)
    )

  plot_size <- max(8, min(14, 6 + num_tips * 0.07))
  ggsave(
    output_image,
    plot = p,
    width = plot_size,
    height = plot_size,
    dpi = 160,
    bg = "white",
    limitsize = FALSE
  )
} else {
  plot_height <- max(8, num_tips * 0.18)
  x_max <- max(ggtree(tree)$data$x, na.rm = TRUE) * 1.35

  p <- ggtree(tree, size = 0.35) +
    geom_tiplab(size = tip_size, color = "#0a160a", offset = 0.002, align = TRUE) +
    geom_nodelab(
      aes(
        subset = !is.na(label) & label != "" & suppressWarnings(as.numeric(label)) >= support_cutoff
      ),
      size = 1.3,
      color = "#1b4332",
      hjust = -0.2
    ) +
    theme_tree2() +
    theme(
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    xlim(0, x_max)

  ggsave(
    output_image,
    plot = p,
    width = 12,
    height = plot_height,
    dpi = 160,
    bg = "white",
    limitsize = FALSE
  )
}

pdf_output <- sub("\\.png$", ".pdf", output_image)
if (layout == "circular") {
  ggsave(pdf_output, plot = p, width = plot_size, height = plot_size, bg = "white", limitsize = FALSE)
} else {
  ggsave(pdf_output, plot = p, width = 12, height = plot_height, bg = "white", limitsize = FALSE)
}

message("Wrote ", output_image, " (", layout, ", ", num_tips, " tips)")
