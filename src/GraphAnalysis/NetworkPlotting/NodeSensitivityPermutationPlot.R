library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)

fd_sens_file = "work/full-drugbank/clustering/metrics/sensitivity.csv"
fd_sens_norm_file= "work/full-drugbank/clustering/metrics/sensitivity_norm.csv"
fd_sa_sens_file = "work/full-drugbank-benchmark/clustering/metrics/sensitivity_sab.csv"
fd_zscore_sens_file = "work/full-drugbank-benchmark/clustering/metrics/sensitivity_zscore.csv"
ma_sens_file = "work/MedAdr/clustering/metrics/sensitivity.csv"
ma_sens_norm_file= "work/MedAdr/clustering/metrics/sensitivity_norm.csv"
ma_sa_sens_file = "work/MedAdr-benchmark/clustering/metrics/sensitivity_sab.csv"
ma_zscore_sens_file = "work/MedAdr-benchmark/clustering/metrics/sensitivity_zscore.csv"
hpo_file= "work/MedAdr/clustering/metrics/hpo_sensitivity.csv"
hpo_norm_file= "work/MedAdr/clustering/metrics/hpo_sensitivity_norm.csv"
hpo_sa_file = "work/MedAdr-benchmark/clustering/metrics/hpo_sensitivity_sab.csv"
hpo_zscore_file= "work/MedAdr-benchmark/clustering/metrics/hpo_sensitivity_zscore.csv"

args <- commandArgs(trailingOnly=TRUE)

fd_sens_file <- args[1]
fd_sens_norm_file <- args[2]
fd_sa_sens_file <- args[3]
fd_zscore_sens_file <- args[4]

fd_sens_df <- read.csv(fd_sens_file, sep = "\t")
fd_sens_norm_df <- read.csv(fd_sens_norm_file, sep = "\t")
fd_sa_sens_df <- read.csv(fd_sa_sens_file, sep = "\t")
fd_zscore_sens_df <- read.csv(fd_zscore_sens_file, sep = "\t")

fd_sens_df$Group = "Drug-nodes: drug-drug Network"
fd_sens_df$Method = "InBEW"
fd_sens_norm_df$Group = "Drug-nodes: drug-drug Network"
fd_sens_norm_df$Method = "Normalised InBEW"
fd_sa_sens_df$Group = "Drug-nodes: drug-drug Network"
fd_sa_sens_df$Method = "Sab"
fd_zscore_sens_df$Group = "Drug-nodes: drug-drug Network"
fd_zscore_sens_df$Method = "Compliment CDF(Zd)"

fd_df <- bind_rows(
  list(
    fd_sens_df,
    fd_sens_norm_df,
    fd_sa_sens_df,
    fd_zscore_sens_df
    ))


ma_sens_file <- args[5]
ma_sens_norm_file <- args[6]
ma_sa_sens_file <- args[7]
ma_zscore_sens_file <- args[8]

ma_sens_df <- read.csv(ma_sens_file, sep="\t")
ma_sens_norm_df <- read.csv(ma_sens_norm_file, sep="\t")
ma_sa_sens_df <- read.csv(ma_sa_sens_file, sep="\t")
ma_zscore_sens_df <- read.csv(ma_zscore_sens_file, sep="\t")

ma_sens_df$Group = "Drug-nodes:drug-ADR Network"
ma_sens_df$Method = "InBEW"
ma_sens_norm_df$Group = "Drug-nodes:drug-ADR Network"
ma_sens_norm_df$Method = "Normalised InBEW"
ma_sa_sens_df$Group = "Drug-nodes:drug-ADR Network"
ma_sa_sens_df$Method = "Sab"
ma_zscore_sens_df$Group = "Drug-nodes:drug-ADR Network"
ma_zscore_sens_df$Method = "Compliment CDF(Zd)"

ma_df <- bind_rows(
  ma_sens_df,
  ma_sens_norm_df,
  ma_sa_sens_df,
  ma_zscore_sens_df
)


hpo_file <- args[9]
hpo_norm_file <- args[10]
hpo_sa_file <- args[11]
hpo_zscore_file <- args[12]


hpo_file_df <- read.csv(hpo_file, sep="\t")
hpo_norm_df <- read.csv(hpo_norm_file, sep="\t")
hpo_sa_df <- read.csv(hpo_sa_file, sep="\t")
hpo_zscore_df <- read.csv(hpo_zscore_file, sep="\t")

hpo_file_df$Group = "HPO-nodes: drug-ADR Network"
hpo_file_df$Method = "InBEW"
hpo_norm_df$Group = "HPO-nodes: drug-ADR Network"
hpo_norm_df$Method = "Normalised InBEW"
hpo_sa_df$Group = "HPO-nodes: drug-ADR Network"
hpo_sa_df$Method = "Sab"
hpo_zscore_df$Group = "HPO-nodes: drug-ADR Network"
hpo_zscore_df$Method = "Compliment CDF(Zd)"

hpo_df <- bind_rows(
  hpo_file_df,
  hpo_norm_df,
  hpo_sa_df,
  hpo_zscore_df
)


output_file <- args[13]


full_df = bind_rows(fd_df, ma_df, hpo_df)

full_df <- full_df[order(full_df$Group, full_df$Sensitivity), ]
sn <- factor(x = 1:nrow(full_df), labels = full_df$Node)
full_df$ord <- sn


g1 <- ggplot(filter(full_df, Group=="Drug-nodes: drug-drug Network"),
             aes(x=ord, y=Sensitivity, color=Method)) +
  geom_jitter(height = 0) +
  ylim(0,1) +
  theme_bw() +
  ggtitle("Drug-node: drug-drug network") +
  xlab("") +
  theme(
    axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='black',),
    legend.position = "none"
  )


g2 <- ggplot(filter(full_df, Group=="Drug-nodes:drug-ADR Network"),
             aes(x=ord, y=Sensitivity, color=Method)) +
  geom_jitter(height = 0) +
  ggtitle("Drug-node: drug-ADR network") +
  ylim(0,1) +
  theme_bw() +
  xlab("") +
  theme(
    axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='black',),
    legend.position = "none"
  )


g3 <- ggplot(filter(full_df, Group=="HPO-nodes: drug-ADR Network"),
             aes(x=ord, y=Sensitivity, color=Method)) +
  geom_jitter(height = 0) +
  ylim(0,1) +
  theme_bw() +
  ggtitle("HPO-node: drug-ADR network") +
  xlab("") +
  theme(
    axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='black',),
    legend.position = "bottom"
  )
  
g <- grid.arrange(g1,g2,g3)


ggsave("work/meta_plots/Sensitivity.png", g, dpi=300, height = 225, width = 170, units = "mm")
