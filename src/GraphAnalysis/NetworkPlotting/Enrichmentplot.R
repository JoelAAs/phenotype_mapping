library(ggplot2)
library(ggrepel)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
enrichment_results <- args[1]
output_file <- args[2]

cluster_enrichment <- read.csv(
    enrichment_results,
    sep = "\t")

if (nrow(cluster_enrichment) != 0) {
  cluster_enrichment <- cluster_enrichment[order(-log10(cluster_enrichment$qvalue), decreasing = T), ]
  cluster_enrichment$show_label <- F
  cluster_enrichment[1:10, "show_label"] <- T


  cluster_enrichment$BG <- sapply(cluster_enrichment$BgRatio, function(x) eval(parse(text = x)))
  cluster_enrichment$sig <- sapply(cluster_enrichment$qvalue, function(x) x < 0.05)

  g <- ggplot(cluster_enrichment, aes(x=BG, y=-log10(qvalue), color=sig)) +
    scale_color_manual(values=c("#c3d3ef", "#276FBF")) +
    theme_bw() +
    geom_point() +
    geom_hline(
      yintercept = -log10(0.05),
      color = "grey40",
      linetype = "dashed") +
    geom_text_repel(
      label = ifelse(cluster_enrichment$show_label & cluster_enrichment$sig, cluster_enrichment$Description, ""),
      max.overlaps = Inf,force=12, point.padding = 2, nudge_y = 3, nudge_x = 0.025, color="black", size=2.5) +
    theme(
      legend.position = "none"
    )


  ggsave(output_file, g , dpi = 300, height = 85, width = 85, units = "mm")
} else {
  # Snakemake wants output
  g <- ggplot() + theme_void()
  ggsave(output_file, g , dpi = 300)
}