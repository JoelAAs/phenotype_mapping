import pandas as pd
import json

rule calculate_term_distances:
    input:
        gene_distances = "work/{project}/gene_distance/all_unique_gene_distance.csv",
        from_to = "work/{project}/gene_interactions/{term_a}_{term_b}.csv.bz2",
        interaction_endpoints= "work/{project}/gene_distance/all_endpoint.json",
        term_a = "{input_location}/{{project}}/{{term_a}}.csv".format(input_location=config['input_location']),
        term_b = "{input_location}/{{project}}/{{term_b}}.csv".format(input_location=config['input_location'])
    output:
        term_distance = "work/{project}/gene_distance/{term_a}_{term_b}.csv"
    run:
        def get_score(df_all, from_genes, depth):
            path_summed = df_all[df_all["from"].isin(from_genes)][df_all["depth"] == depth].paths.sum()
            tot_paths = 0
            for gene in from_genes:
                tot_paths += endpoints[gene][str(depth)]["possible_paths"]
            return path_summed / tot_paths

        df_distance = pd.read_csv(input.gene_distances, sep="\t")
        df_from_to = pd.read_csv(input.from_to, sep="\t")
        df_rev = df_from_to.rename(columns={"gene_a": "gene_b", "gene_b": "gene_a"})
        df_from_to = pd.concat([df_from_to, df_rev])
        df_from_to = df_from_to[~df_from_to.duplicated()]
        df_from_to = df_from_to.rename(columns={"gene_a": "from", "gene_b": "to"})
        df_all = df_distance.merge(df_from_to,how="right",on=["from", "to"])
        endpoints = json.loads([l for l in open(input.interaction_endpoints,"r")][0])

        genes_a = list(pd.read_csv(input.term_a).gene)
        genes_b = list(pd.read_csv(input.term_b).gene)
        with open(output.term_distance, "w") as w:
            w.write("from_term\tdepth\tprobability\n")
            for i in range(1, config["max_depth"]):
                prob_a = get_score(df_all, genes_a, i)
                w.write(f"{wildcards.term_a}\t{i}\t{prob_a}\n")
                prob_b = get_score(df_all, genes_b, i)
                w.write(f"{wildcards.term_b}\t{i}\t{prob_b}\n")
