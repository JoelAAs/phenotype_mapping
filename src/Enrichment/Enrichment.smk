
rule KEGGEnchriment:
    input:
        "work/edgelists/clustering/candidate_gene/inferred/{joinedname}_{cutoff}_{method}_{cluster}_MCODE_entrez.csv"
    output:
        "work/edgelists/clustering/candidate_gene/enrichment/{joinedname}_{cutoff}_{method}_{cluster}_KEGG.csv"
    shell:
        """
        Rscript src/Enrichment/Enrichment.R {input} {output}
        """