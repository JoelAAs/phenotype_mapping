import os
import glob
import pandas as pd
import networkx as nx
import random

drugs = [os.path.basename(f).replace(".csv", "")
         for f in  glob.glob(f"input/full-drugbank/*.csv")]
permuts = 100


rule GenerateRandomSets:
    params:
        limit = 0.7
    input:
        ppi_file = "data/9606.protein.links.v11.5.txt",
        term = "input/full-drugbank/{drug}.csv"
    output:
        permuts = expand("work/full-drugbank/permut/{{drug}}_set_{n}.csv", n = range(permuts))
    run:
        edge_list_df = pd.read_csv(input.ppi_file, sep=" ")
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

        for i in range(permuts):
            output_file = f"work/full-drugbank/permut/{wildcards.drug}_set_{i}.csv"
            with open(output_file, "w") as w:
                w.write("gene\n")
                for gene in genes:
                    selected_gene = random.sample(
                        degree_dict[node_degree_dict[gene]],1)[0]
                    w.write(f"{selected_gene}\n")


rule unique_genes_perumt:
    input:
        permut_sets = expand("work/full-drugbank/permut/{drug}_set_{n}.csv", n = range(permuts), drug=drugs)
    output:
        unique_genes = "work/full-drugbank/permut/unique_genes.csv"
    run:
        unique_genes = set()
        for genefile in input.permut_sets:
            with open(genefile, "r") as f:
                c_set = {gene.strip() for gene in f.readlines()[1:]}
                unique_genes.update(c_set)

        with open(output.unique_genes, "w") as w:
            w.write("gene\n")
            for gene in unique_genes:
                w.write(f"{gene}\n")
