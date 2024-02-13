import math

import numpy as np
import pandas as pd
import networkx as nx

rule get_zscores:
    params:
        ppi_file = "data/9606.protein.links_above_700.v11.5.txt",
        limit = 0.7
    input:
        gene_probabilities = "data/full_drugbank_gene_probabilites/{group}_gene_probabilites.csv",
        set_permutations = "work/full-drugbank-benchmark/candidate_genes/term_probabilities_groups/{group}_mean_var.csv"
    output:
        z_scores = "work/full-drugbank-benchmark/group-z-scores/{group}_z_score.csv"
    run:
        n_permuts = 100 # TODO: config later
        set_permut_dict = dict()
        with open(input.set_permutations, "r") as w:
            for l in w.readlines()[1:]:
                gene, mean, var = l.strip().split("\t")
                set_permut_dict[gene] = [float(mean), float(var)]

        # PPI Degree conf
        edge_list_df = pd.read_csv(params.ppi_file,sep=" ")
        edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
        edge_list_df = edge_list_df[edge_list_df["combined_score"] >= params.limit]
        G = nx.from_pandas_edgelist(
            edge_list_df,
            "protein1",
            "protein2"
        )

        degree_dict = {}
        for node, c_degree in G.degree():
            if c_degree in degree_dict:
                degree_dict[c_degree].append(node)
            else:
                degree_dict[c_degree] = [node]
        node_degree_dict = {node: degree for node, degree in G.degree()}

        mean_var_dict = dict()
        def get_sample_mean_var(degree):
            if degree not in mean_var_dict:
                genes = degree_dict[degree]
                i = 1
                while len(genes) < 4:
                    if degree-i in degree_dict:
                        genes += degree_dict[degree-i]
                    i += 1
                means = np.zeros(len(genes))
                variances = np.zeros(len(genes))
                for i, gene in enumerate(genes):
                    if gene in set_permut_dict:
                        c_m, c_v = set_permut_dict[gene]
                        means[i] = math.log10(c_m)
                        variances[i] = math.log10(c_v)

                degree_mean = sum(means)/len(means)
                degree_variance = means.var()
                mean_var_dict[degree] = [degree_mean, degree_variance]
            return mean_var_dict[degree]

        with open(input.gene_probabilities, "r") as f:
            with open(output.z_scores, "w") as w:
                w.write("Gene\ty_probability\tz_score\tdegree\n")
                lines = f.readlines()[1:]
                genes = [l.strip().split("\t")[0] for l in lines]
                probs = [float(l.strip().split("\t")[1]) for l in lines]
                for gene, prob in zip(genes, probs):
                    if gene not in node_degree_dict:
                        w.write(f"{gene}\t{prob}\t-\t-\n")
                    else:
                        mu, var = get_sample_mean_var(node_degree_dict[gene])
                        z_score = (math.log10(prob) - mu)/math.sqrt(var)
                        w.write(f"{gene}\t{prob}\t{z_score}\t{node_degree_dict[gene]}\n")

                not_reached = [gene for gene in node_degree_dict.keys() if gene not in genes]
                for gene_left in not_reached:
                    mu, var = get_sample_mean_var(node_degree_dict[gene_left])
                    z_score = (0 - mu)/math.sqrt(var)
                    w.write(f"{gene_left}\t{0}\t{z_score}\t{node_degree_dict[gene_left]}\n")


rule get_top_values:
    params:
        n_top = 200
    input:
        z_scores = "work/full-drugbank-benchmark/group-z-scores/{group}_z_score.csv",
        unique_input = "data/full_drugbank_gene_probabilites/{group}_unique_genes.csv"
    output:
        top_z_scores = "work/full-drugbank-benchmark/group-z-scores/{group}_top.csv"
    run:
        z_scores = pd.read_csv(input.z_scores, sep="\t")
        print(z_scores)
        with open(input.unique_input, "r") as f:
            input_genes = [g.strip() for g in f.readlines()[1:]]
        z_scores = z_scores[~z_scores["Gene"].isin(input_genes)]
        z_scores_top = z_scores.sort_values(by="z_score", ascending=False)[:params.n_top]
        z_scores_top.to_csv(output.top_z_scores, index=False, sep = "\t")