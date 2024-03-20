library(ggplot2)
groups <- c("Statin", "Antidepressant", "NSAID")

prob_dfs = list()
i = 1
for (group in groups) {
  df_prob <- read.csv(
    paste0("work/full-drugbank-benchmark/group-quant/", group, "_quant.csv"),
    sep="\t")
  df_diamond <- read.csv(
    paste0("work/full-drugbank-benchmark/benchmarking/DIAMOnD_", group, ".csv"),
    sep="\t")
  df_diamond["Selection"] <- "DIAMOnD"
  df_diamond <- df_diamond[, c("X.rank" ,"DIAMOnD_node", "Selection")]
  df_prob <- merge(df_prob, df_diamond, by.x="Gene", by.y="DIAMOnD_node", all.x = T)
  df_prob["Group"] <- group
  prob_dfs[[i]] <- df_prob
  i = i + 1
}

all_prob <- do.call(rbind, prob_dfs)
all_prob[is.na(all_prob$Selection), "Selection"] = "CDF"

g1 <-  ggplot(all_prob,
       aes(y=quant, x=X.rank)
       ) +
  geom_bin2d() +
  scale_fill_gradient(low="lightblue", high="red") +
 theme_bw() +
  xlab("DIAMOnD Rank") +
  ylab("CDF(P(y|C))") +
#  ggtitle("CDF of P(y|C) vs DIAMOnD rank") +
  facet_grid(. ~Group )

ggsave("paper/manual_figures/diamond_vs_CDF.png", g1, dpi=300, height = 60, width = 170, units = "mm")
