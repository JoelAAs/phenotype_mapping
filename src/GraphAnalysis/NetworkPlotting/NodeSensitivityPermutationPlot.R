library(ggplot2)
library(ggrepel)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
fd_sens_file <- args[1]


ma_sens_file <- args[2]


hpo_file <- args[3]


output_file <- args[4]

fd_sens_df <- read.csv(fd_sens_file, sep = "\t")
fd_sens_df$Group = "Drug-nodes in Drug Network"
ma_sens_df <- read.csv(ma_sens_file, sep = "\t")
ma_sens_df$Group = "Drug-nodes in Drug-ADR Network"
hpo_df <- read.csv(hpo_file, sep = "\t")
hpo_df$Group = "ADR-nodes in Drug-ADR Network"

full_df = rbind(hpo_df, ma_sens_df, fd_sens_df)

full_df <- full_df[order(full_df$Group, full_df$Sensitivity), ]
sn <- factor(x = 1:nrow(full_df), labels = full_df$Node)
full_df$ord <- sn

g <- ggplot(full_df, aes(x=ord, y=Sensitivity, fill=Group)) +
  geom_col(width = 0.85) +
  labs() +
  labs(shape='From term genesets',
       colour="Within gensets (%)") +
  xlab("Node within cluster") +
  ylab("Sensitivity") +
  facet_grid(~Group, scales="free_x") +
  theme(
    axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='black',),
    legend.position = "none"
    )


ggsave(output_file, g, dpi=300)