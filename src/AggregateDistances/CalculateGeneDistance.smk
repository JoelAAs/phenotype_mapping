import pandas as pd
import bz2
import os
import json

#### Input Functions
def get_all_genes_at_depth(wc):

    all_needed_files_exists = all(
        map(
            os.path.exists,
        [
            f"work/{wc.project}/gene_path",
            f"work/{wc.project}/gene_interactions/unique_genes.csv"]))
    if not all_needed_files_exists:
        re_eval = checkpoints.get_possible_paths_from_genes.get(**wc)
        return re_eval
    else:
        df_distance_pairs = pd.read_csv(f"work/{wc.project}/gene_interactions/unique_genes.csv", sep = "\t")
        expected_distance = []
        for i, gene in df_distance_pairs.itertuples():
            for depth_expected in range(1, config["max_depth"] +1):
                expected_distance.append(
                    f"work/{wc.project}/gene_path/{gene}_at_{depth_expected}_sum.csv.bz2"
                )

        return expected_distance


#### Rules
rule get_list_of_all_gene_depths:
    input:
        genes = lambda wildcards: get_all_genes_at_depth(wildcards)
    output:
        expected_output = "work/{project}/checkpoint/all_genes_at_depth.txt"
    run:
        with open(output.expected_output, "w") as w:
            for gene in input.genes:
                w.write(f"{gene}\n")


rule create_depth_target_json:
    input:
        genes = "work/{project}/gene_interactions/unique_genes.csv",
        checkpoint = "work/{project}/checkpoint/all_genes_at_depth.txt"
    output:
        interaction_endpoints = "work/{project}/gene_distance/all_endpoint.json"
    run:
        output_dict = dict()
        with open(input.genes, "r") as f:
            lines = [l.strip() for l in f]
        genes = lines[1:]

        for gene in genes:
            output_dict[gene] = dict()
            for depth in range(1, config["max_depth"] + 1):
                output_dict[gene][depth] = dict()
                with bz2.open(f"work/{wildcards.project}/gene_path/{gene}_at_{depth}_sum.csv.bz2","r") as f:
                    endpoint_lines = f.readlines()
                    endpoint_lines = [l.decode("utf-8").strip() for l in endpoint_lines]
                    number_paths = 0
                    for gene_path_line in endpoint_lines[1:]:
                        gene_path_line = gene_path_line.split()
                        endpoint_gene = gene_path_line[0]
                        number_of_hits = int(gene_path_line[1])
                        number_paths += number_of_hits
                        output_dict[gene][depth][endpoint_gene] = number_of_hits
                    output_dict[gene][depth]["possible_paths"] = number_paths

        with open(output.interaction_endpoints, "w") as w:
            w.write(json.dumps(output_dict))


rule calculate_gene_distance:
    input:
        all_unique_pairs = "work/{project}/gene_interactions/unique.csv",
        interaction_endpoints = "work/{project}/gene_distance/all_endpoint.json"
    output:
        distance = "work/{project}/gene_distance/all_unique_gene_distance.csv"
    run:
        def _get_depth(gene, target, endpoints):
            for i in range(1, config["max_depth"] + 1):
                i = str(i)
                if target in endpoints[gene][i]:
                    number_of_hits = endpoints[gene][i][target]
                    return number_of_hits, i , endpoints[gene][i]["possible_paths"]

            return 0, -1

        endpoints = json.loads([l for l in open(input.interaction_endpoints, "r")][0])

        df_unique = pd.read_csv(input.all_unique_pairs, sep = "\t")
        tot_rows = len(df_unique)
        with open(output.distance, "w") as w:
            w.write("from\tto\tdepth\tpaths\tpossible_path\n")
            for i, gene_a, gene_b in df_unique.itertuples():
                hit_a, depth_a, possble_paths_a = _get_depth(gene_a, gene_b, endpoints)
                w.write(f"{gene_a}\t{gene_b}\t{depth_a}\t{hit_a}\t{possble_paths_a}\n")
                hit_b, depth_b, possble_paths_b = _get_depth(gene_b, gene_a, endpoints)
                w.write(f"{gene_b}\t{gene_a}\t{depth_b}\t{hit_b}\t{possble_paths_b}\n")
                print(f"{i}/{tot_rows} rows: {round(i/tot_rows, 3)*100} % done", end="\r")



