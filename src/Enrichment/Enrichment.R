library(clusterProfiler)

### Main
args <- commandArgs(trailingOnly=TRUE)
genes_entrez_file <- args[1]
output_file <- args[2]

genes <- read.csv(genes_entrez_file, sep="\t")
genes <- genes$entrez

pathEnrich <- enrichKEGG(
    gene = genes,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH")

write.csv(pathEnrich@result, output_file, sep="\t")
