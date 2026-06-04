#!/usr/bin/env Rscript
# Render phylogenetic trees for FungiFlow (rectangular or fan/circular).
# Usage: Rscript plot_tree.R <treefile> <output.png> <layout> [tip_meta.csv]
#   tip_meta.csv columns: cluster, genus, species, tip_label

suppressPackageStartupMessages({
  library(ggtree)
  library(ggplot2)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_tree.R <treefile> <output.png> [layout] [tip_meta.csv]")
}

tree_file <- args[1]
output_image <- args[2]
layout <- if (length(args) >= 3) tolower(args[3]) else "rectangular"
meta_file <- if (length(args) >= 4) args[4] else ""

if (!file.exists(tree_file)) {
  stop("Tree file not found: ", tree_file)
}

tree <- read.tree(tree_file)
num_tips <- length(tree$tip.label)
tree$node.label <- as.character(tree$node.label)

tip_size <- if (num_tips > 80) 0.9 else if (num_tips > 40) 1.1 else 1.5
support_cutoff <- 70

fungi_colors <- c(
  "#1b4332", "#2d6a4f", "#40916c", "#52b788", "#74c69d", "#95d5b2",
  "#bc6c25", "#9b2226", "#4a6741", "#3d405b", "#81b29a", "#e07a5f",
  "#6a4c93", "#1982c4", "#8ac926", "#ff595e"
)

load_tip_meta <- function(tree) {
  # ggtree %<+% joins on a single column named "label" — do not add a second one.
  tip_df <- data.frame(
    label = tree$tip.label,
    tip_label = tree$tip.label,
    genus = "Unknown",
    species = "",
    stringsAsFactors = FALSE
  )

  if (!nzchar(meta_file) || !file.exists(meta_file)) {
    return(tip_df)
  }

  meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cluster" %in% names(meta)) {
    return(tip_df)
  }

  meta$cluster <- as.character(meta$cluster)

  match_cluster <- function(cid) {
    if (cid %in% meta$cluster) {
      return(cid)
    }
    hits <- meta$cluster[
      meta$cluster == cid |
        endsWith(cid, meta$cluster) |
        endsWith(meta$cluster, cid)
    ]
    if (length(hits) > 0) {
      return(hits[1])
    }
    NA_character_
  }

  for (i in seq_len(nrow(tip_df))) {
    cid <- tip_df$label[i]
    hit <- match_cluster(cid)
    if (is.na(hit)) {
      next
    }
    row <- meta[meta$cluster == hit, , drop = FALSE][1, , drop = FALSE]
    if ("tip_label" %in% names(row) && nzchar(row$tip_label[1])) {
      tip_df$tip_label[i] <- row$tip_label[1]
    }
    if ("genus" %in% names(row) && nzchar(row$genus[1])) {
      tip_df$genus[i] <- row$genus[1]
    }
    if ("species" %in% names(row) && nzchar(row$species[1])) {
      tip_df$species[i] <- row$species[1]
    }
  }

  tip_df$genus <- ifelse(is.na(tip_df$genus) | tip_df$genus == "", "Unknown", tip_df$genus)
  tip_df
}

tip_df <- load_tip_meta(tree)
n_genera <- length(unique(tip_df$genus))
pal <- rep(fungi_colors, length.out = max(n_genera, 1))

p_base <- function() {
  ggtree(tree, aes(color = genus), size = 0.35) %<+% tip_df +
    scale_color_manual(
      name = "Genus",
      values = pal,
      na.value = "#7A7A7A"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, linewidth = 1)))
}

# Fan: open_tree() leaves an empty sector (not a full 360° ring).
fan_open_angle <- 120

if (layout == "circular") {
  p <- ggtree(tree, layout = "circular", aes(color = genus), size = 0.35) %<+% tip_df +
    scale_color_manual(
      name = "Genus",
      values = pal,
      na.value = "#7A7A7A"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, linewidth = 1)))
  p <- open_tree(p, angle = fan_open_angle)
  p <- p +
    geom_tiplab2(
      aes(label = tip_label, color = genus, angle = angle),
      size = tip_size,
      align = TRUE,
      offset = 0.003,
      linewidth = 0.15,
      linecolor = "grey70"
    ) +
    geom_nodelab(
      aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= support_cutoff),
      size = 1.1,
      color = "#1b4332"
    ) +
    theme(
      plot.margin = margin(10, 10, 10, 10),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    )

  plot_size <- max(9, min(16, 7 + num_tips * 0.08))
  ggsave(
    output_image,
    plot = p,
    width = plot_size,
    height = plot_size,
    dpi = 170,
    bg = "white",
    limitsize = FALSE
  )
} else {
  plot_height <- max(9, num_tips * 0.22)
  x_max <- max(ggtree(tree)$data$x, na.rm = TRUE) * 1.55

  p <- p_base() +
    geom_tiplab(
      aes(label = tip_label, color = genus),
      size = tip_size,
      offset = 0.004,
      align = TRUE,
      linecolor = "grey70",
      linewidth = 0.2
    ) +
    geom_nodelab(
      aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= support_cutoff),
      size = 1.2,
      color = "#1b4332",
      hjust = -0.15
    ) +
    theme_tree2() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    ) +
    xlim(0, x_max)

  ggsave(
    output_image,
    plot = p,
    width = 14,
    height = plot_height,
    dpi = 170,
    bg = "white",
    limitsize = FALSE
  )
}

pdf_output <- sub("\\.png$", ".pdf", output_image)
if (layout == "circular") {
  ggsave(pdf_output, plot = p, width = plot_size, height = plot_size, bg = "white", limitsize = FALSE)
} else {
  ggsave(pdf_output, plot = p, width = 14, height = plot_height, bg = "white", limitsize = FALSE)
}

fan_note <- if (layout == "circular") {
  paste0("fan open=", fan_open_angle)
} else {
  layout
}
message("Wrote ", output_image, " (", fan_note, ", ", num_tips, " tips)")
