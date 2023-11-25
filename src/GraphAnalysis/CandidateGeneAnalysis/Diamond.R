library(MODifieR)

read_ppi <- function(ppi_file, cutoff) {
        ppi_df <- read.delim(
            ppi_file,
            stringsAsFactors=FALSE,
            sep=" "
        )
        ppi_df <- ppi_df[ppi_df$combined_score > cutoff*1000, ]
        return(ppi_df)
}

read_genes <- function(hit_genes, top_percent) {
    top_df <- read.delim(
        hit_genes,
        stringsAsFactors=FALSE
    )
    top_df$gene <- top_df$entrez
    top_df$pvalue <- 0
    top_df <- top_df[, c("gene", "pvalue")]
    hit_input = create_custom_microarray_input_object(
        diff_genes=top_df
    )
    return(hit_input)
}


### Main
args <- commandArgs(trailingOnly=TRUE)
hit_file <- args[1]
output_file <- args[2]
ppi_file <- args[3]
ppi_cutoff <- as.numeric(args[4])

ppi_df <- read_ppi(ppi_file, ppi_cutoff)
hit_input <- read_genes(hit_file, top_percent)

if (nrow(hit_input$diff_genes) > 1) {
    module_diamond <- diamond(
        hit_input,
        ppi_df,
        include_seed=T
        )

    diamond_gene_df <- as.data.frame(module_diamond$module_genes, col.names=c("genes"))
    colnames(diamond_gene_df) <- c("gene")
    write.table(diamond_gene_df, output_file, row.names = FALSE, sep = "\t")
} else {
    write("gene", output_file)
}