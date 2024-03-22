import os
import glob
import pandas as pd
import networkx as nx
import random

terms = [os.path.basename(f).replace(".csv", "")
         for f in  glob.glob(f"input/{config['permutation_folder']}/*.csv")]
n_permutations = config["permutation_N"]
input_folder = config["project_name"]

rule GenerateRandomSets:
    params:
        limit = config["combined_score_threshold"],
        ppi_network_file= config["ppi_file"]
    input:
        term = "input/{input_folder}/{{term}}.csv".format(input_folder=input_folder)
    output:
        permuts = expand("input/{permutation_folder}/{{term}}_set_{n}.csv", permutation_folder=config['permutation_folder'], n = range(n_permutations))
    run:
        edge_list_df = pd.read_csv(params.ppi_network_file, sep=" ")
        edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
        edge_list_df = edge_list_df[edge_list_df["combined_score"] >= params.limit]
        G = nx.from_pandas_edgelist(
            edge_list_df,
            "protein1",
            "protein2",
            edge_attr=True
        )
        node_degree_dict = {node: degree for node, degree in G.degree()}
        degree_dict = {}
        for node, c_degree in G.degree():
            if c_degree in degree_dict:
                degree_dict[c_degree].append(node)
            else:
                degree_dict[c_degree] = [node]

        with open(input.term, "r") as f:
            genes = {gene.strip() for gene in f.readlines()[1:]}

        genes = genes & G.nodes()

        for i in range(n_permutations):
            output_file = f"input/{config['permutation_folder']}/{wildcards.term}_set_{i}.csv"
            with open(output_file, "w") as w:
                w.write("gene\n")
                for gene in genes:
                    selected_gene = random.sample(
                        degree_dict[node_degree_dict[gene]],1)[0]
                    w.write(f"{selected_gene}\n")


rule GeneratePermutationYaml:
    output:
        yaml = "input/{permutation_folder}.yaml"
    run:
        with open(output.yaml, "w") as w:
            w.write(f"project_name: {config['permutation_folder']}\n")
            w.write(f"ppi_file: \"{config['ppi_file']}\"\n")
            w.write(f"combined_score_threshold: {config['combined_score_threshold']}\n")
            w.write(f"max_depth: {config['max_depth']}\n")
            w.write(f"n_path_batches: {config['n_path_batches']}\n")
            w.write(f"n_permutation_combinations: {config['n_permutation_combinations']}\n")