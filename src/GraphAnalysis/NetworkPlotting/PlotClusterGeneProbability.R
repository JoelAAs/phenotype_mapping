library(ggplot2)
library(ggrepel)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
probabilities <- args[1]
output_file <- args[2]

prob_df <- read.csv(
  probabilities,
  sep = "\t")

prob_df_threshold <- mean(prob_df$y_probability) + 2*sd(prob_df$y_probability)
prob_df$shape <- sapply(prob_df$count, function(x) x == 0)
prob_df$show_label <- sapply(prob_df$y_probability, function(x) x> prob_df_threshold)
prob_df$centrality_norm <- prob_df$centrality /max(prob_df$centrality)

prob_df$presence <- prob_df$count/prob_df$n_terms
g <- ggplot(prob_df, aes(y=centrality_norm, x=log(y_probability), color=presence, shape=shape)) +
  scale_shape_discrete(labels=c('No', 'Yes'),) +
  labs(shape='From term genesets',
       colour="Within gensets (%)") +
  xlab("Log probability of gene within cluster") +
  ylab("Normalized degree centrality of gene within PPI") +
  theme_bw() +
  geom_vline(
    xintercept = log(prob_df_threshold),
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
  theme(
    axis.text=element_text(size=8),
    axis.title=element_text(size=16),
    plot.title=element_text(size=20),
    strip.text = element_text(size=12)
  )

ggsave(output_file, g, dpi=300)