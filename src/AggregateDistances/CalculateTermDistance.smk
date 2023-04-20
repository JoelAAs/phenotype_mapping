import pandas as pd

rule calculate_term_distances:
    input:
        gene_distances = "work/{project}/gene_distance/all_unique_gene_distance.csv",
        from_to = "work/{project}/gene_interactions/{term_a}_{term_b}.csv.bz2"
    output:
        term_distance = "work/{project}/gene_distance/{term_a}_{term_b}.csv.bz2"
    run:
        df_distance = pd.read_csv(input.gene_distances, sep="\t")
        df_from_to = pd.read_csv(input.from_to, sep="\t")

        df_all_first = df_distance.merge(df_from_to, how="right", on=["gene_a", "gene_b"])
        df_all_flip = df_all_first[(df_all_first.depth != df_all_first.depth)]  #get all without hits
        df_all_flip = df_all_flip[[
            "gene_a",
            "gene_b"
        ]]
        df_all_flip.columns = [
            "gene_b",
            "gene_a"
        ] # Flip
        df_all_flip = df_distance.merge(df_all_flip,how="right",on=["gene_a", "gene_b"])
        df_all = pd.concat([df_all_first, df_all_flip])
        df_all = df_all[(df_all.depth == df_all.depth)]
        df_all.to_csv(output.term_distance, sep="\t", index = False)