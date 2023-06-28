import pandas as pd

rule MCODE:
    params:
        ppi_cutoff=0.7,
        top_percent=0.05
    input:
        cluster_genes = "work/edgelists/clustering/candidate_gene/{joinedname}_{cutoff}_{method}_{cluster}.csv",
        ppi_cluster = "data/9606.protein.links.v11.5.txt"
    output:
        inferred_genes = "work/edgelists/clustering/candidate_gene/inferred/{joinedname}_{cutoff}_{method}_{cluster}_MCODE.csv"
    singularity:
        "envs/modifier.simg"
    shell:
        """
        Rscript src/Enrichment/MCODE.R {input.cluster_genes} {output.inferred_genes} \
            {input.ppi_cluster} {params.ppi_cutoff} {params.top_percent} 
        """


rule TranslateToGeneName:
    input:
        inferred_genes = "work/edgelists/clustering/candidate_gene/inferred/{joinedname}_{cutoff}_{method}_{cluster}_MCODE.csv",
        stringID = "data/stringdb/9606.protein.info.v11.5.txt"
    output:
        translate_genes = "work/edgelists/clustering/candidate_gene/inferred/{joinedname}_{cutoff}_{method}_{cluster}_MCODE_geneNames.csv"
    run:
        ig_df = pd.read_csv(input.inferred_genes)
        ig_df = ig_df.rename({"genes":"string_protein_id"}, axis =1)
        string_df = pd.read_csv(input.stringID, sep = "\t")
        ig_df = ig_df.merge(string_df, on = "string_protein_id", how ="left")
        ig_df[["preferred_name"]].to_csv(output.translate_genes, index=False)


rule GetEntrezGene:
    input:
        translate_genes = "work/edgelists/clustering/candidate_gene/inferred/{joinedname}_{cutoff}_{method}_{cluster}_MCODE_geneNames.csv",
        entrez = "data/ncbi/entrez.csv"
    output:
        entrez_ids = "work/edgelists/clustering/candidate_gene/inferred/{joinedname}_{cutoff}_{method}_{cluster}_MCODE_entrez.csv"
    run:
        trans_df = pd.read_csv(input.translate_genes)
        trans_df = trans_df.rename({"preferred_name": "gene_name"}, axis=1)
        en_df = pd.read_csv(input.entrez, sep="\t")
        trans_df = trans_df.merge(en_df, on="gene_name", how="left")
        trans_df = trans_df[~trans_df.entrez.isna()]
        trans_df.entrez = [int(eid) for eid in trans_df.entrez]
        trans_df.to_csv(output.entrez_ids, sep="\t", index=False)




