library(dplyr)
library(ggplot2)
all_probs = list()
i = 1
for (project in c("full-drugbank", "MedAdr") ){
  for (cluster in 0:2) {
    file_name = paste0("work/", project, "/candidate_genes/annotated_3/annotated_", cluster, ".csv")
    prob_df <- read.csv(
      file_name,
      sep = "\t")
    prob_df$cluster <- paste0("Cluster: ",cluster)
    if (project == "full-drugbank"){
      prob_df$project <- "Drug-Drug network"
    } else {
      prob_df$project <- "Drug-ADR network"
    }
    prob_df$cutoff <- mean(prob_df$y_probability) + 2*sd(prob_df$y_probability)
    prob_df$centrality_norm <- prob_df$centrality /max(prob_df$centrality)
    all_probs[[i]] <- prob_df
    i = i + 1
  }
}


prob_df <- bind_rows(all_probs)
prob_df$shape <- sapply(prob_df$count, function(x) x == 0)

prob_df$presence <- prob_df$count/prob_df$n_terms
g <- ggplot(prob_df, aes(y=centrality_norm, x=log(y_probability), color=presence, shape=shape)) +
  scale_shape_discrete(labels=c('No', 'Yes'),) +
  labs(shape='From term genesets',
       colour="Within gensets (%)") +
  xlab("Log probability of gene within cluster") +
  ylab("Normalized degree centrality of gene within PPI") +
  theme_bw() +
  geom_vline(
    aes(xintercept = log(cutoff)),
    color = "grey40",
    linetype = "dashed") +
  geom_hline(
    yintercept = 0.3,
    color = "grey40",
    linetype = "dashed") +
  geom_point() +
  scale_colour_gradientn(
    colours = c("#e70040","#785BFF","#609faa", "#163B46", "black"),
    values = c(1.0,0.66,0.33,0.01,0)) +
  facet_grid(project ~ cluster) +
  theme(
    axis.text=element_text(size=8),
    axis.title=element_text(size=16),
    plot.title=element_text(size=20),
    strip.text = element_text(size=12)
  )

ggsave("paper/manual_figures/candidate_dist.png", g, dpi=300)
