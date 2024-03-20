library(ggplot2)
library(ggrepel)
library(dplyr)


all_probs = list()
i = 1
for (project in c("full-drugbank", "MedAdr") ){
  for (cluster in 0:2) {
    
    file_name = paste0("work/", project, "/candidate_genes/enrichment_3/top/enrichment_", cluster, ".csv")
    cluster_enrichment <- read.csv(
      file_name,
      sep = "\t")
    if (cluster == 0) {
      label = "NSAID"
    } else if (cluster == 1) {
      label = "Statin"
    } else {
      label = "Antidepressant"
    }
    cluster_enrichment$cluster <- paste0("Cluster: ", label)
    if (project == "full-drugbank"){
      cluster_enrichment$project <- "Drug-Drug network"
    } else {
      cluster_enrichment$project <- "Drug-ADR network"
    }
    cluster_enrichment <- cluster_enrichment[order(-log10(cluster_enrichment$qvalue), decreasing = T), ]
    cluster_enrichment$show_label <- F
    cluster_enrichment[1:5, "show_label"] <- T
    all_probs[[i]] <- cluster_enrichment
    i = i + 1
  }
}

cluster_enrichment <- bind_rows(all_probs)

cluster_enrichment$BG <- sapply(cluster_enrichment$BgRatio, function(x) eval(parse(text = x)))
cluster_enrichment$sig <- sapply(cluster_enrichment$qvalue, function(x) x < 0.05)

g <- ggplot(cluster_enrichment, aes(x=Count, y=-log10(qvalue), color=sig)) +
  scale_color_manual(values=c("#c3d3ef", "#276FBF")) +
  theme_bw() +
  geom_point() +
  xlab("Genes related to KEGG-pathway") +
  geom_hline(
    yintercept = -log10(0.05),
    color = "grey40",
    linetype = "dashed") +
  geom_text_repel(
    label = ifelse(cluster_enrichment$show_label & cluster_enrichment$sig, cluster_enrichment$Description, ""),
    max.overlaps = 4,force=60, point.padding = 2, nudge_y = 5, nudge_x = 0.025,color="black") +
  facet_grid(project ~ cluster) +
  theme(
    axis.text=element_text(size=8),
    axis.title=element_text(size=16),
    plot.title=element_text(size=20),
    strip.text = element_text(size=12),
    legend.position = "none"
  )


ggsave("paper/manual_figures/Enrichment_all.png", g, width=12, height = 12, dpi = 300)
