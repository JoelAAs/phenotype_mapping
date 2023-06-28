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
    hit_df <- read.delim(
        hit_genes,
        stringsAsFactors=FALSE
    )
    hit_df <- hit_df[order(hit_df$summed_score, decreasing=T), ]
    top_df <- hit_df[1:floor(nrow(hit_df)*top_percent), ]
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
top_percent <- as.numeric(args[5])

ppi_df <- read_ppi(ppi_file, ppi_cutoff)
hit_input <- read_genes(hit_file, top_percent)

module_mcode <- mcode(
    MODifieR_input = hit_input,
    ppi_network = ppi_df,
    haircut = FALSE,
    fluff = False,
    loops = TRUE)

mcode_gene_df <- as.data.frame(module_mcode$module_genes, col.names=c("genes"))
colnames(mcode_gene_df) <- c("genes")
write.table(mcode_gene_df, output_file, row.names = FALSE)