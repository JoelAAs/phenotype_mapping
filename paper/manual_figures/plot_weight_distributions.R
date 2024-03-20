library(ggplot2)
library(dplyr)

fd_inbew <- read.csv("work/full-drugbank/term_to_term_probability_matrix.csv", sep="\t")
fd_inbew$Method <- "InBEW"
fd_norm_inbew <- read.csv("work/full-drugbank/term_to_term_probability_matrix_norm.csv", sep="\t")
fd_norm_inbew$Method <- "Normalised InBEW"
fd_sa <- read.csv("work/full-drugbank-benchmark/term_to_term_probability_matrix_sab.csv", sep="\t")
fd_sa$Method <- "Sab"
fd_za <- read.csv("work/full-drugbank-benchmark/term_to_term_probability_matrix_zscore.csv", sep="\t")
fd_za$Method <- "Compliment CDF(Zd)"

fd_full <- bind_rows(
  fd_inbew,
  fd_norm_inbew,
  fd_sa,
  fd_za
)
fd_full$Network <- "Drug weights"

hpo_inbew <- read.csv("work/HPO-pruned/term_to_term_probability_matrix.csv", sep = "\t")
hpo_inbew$Method <- "InBEW"
hpo_norm_inbew <- read.csv("work/HPO-pruned/term_to_term_probability_matrix_norm.csv", sep="\t")
hpo_norm_inbew$Method <- "Normalised InBEW"
hpo_sa <- read.csv("work/full-drugbank-benchmark/term_to_term_probability_matrix_sab.csv", sep="\t")
hpo_sa$Method <- "Sab"
hpo_za <- read.csv("work/full-drugbank-benchmark/term_to_term_probability_matrix_zscore.csv", sep="\t")
hpo_za$Method <- "Compliment CDF(Zd)"

hpo_full <- bind_rows(
  hpo_inbew,
  hpo_norm_inbew,
  hpo_sa,
  hpo_za
)

hpo_full$Network <- "ADR weights"
all_weights <- rbind(hpo_full, fd_full)

g <- ggplot(all_weights, aes(x=probability, fill=Method)) + 
  geom_histogram() +
  facet_grid(Network ~ Method, scales = "free") +
  ggtitle("Edge-weight distribution per method and network") +
  xlab("Edge-weight") +
  theme_bw() +
  theme(
    axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='black',),
    legend.position = "bottom"
  )

ggsave("paper/manual_figures/Edge_distribution.png", g, dpi=300)



