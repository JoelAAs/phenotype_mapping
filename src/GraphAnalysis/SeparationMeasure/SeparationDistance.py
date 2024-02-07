import os.path

import numpy as np
import networkx as nx
import random
import pandas as pd
import itertools


class CalculateSeparation:
    def __init__(self, ppi_file):
        self.ppi_graph = self.graph_setup(ppi_file)
        self.degree_dict = self.get_degree_dict()
        self.node_degree_dict = {node: degree for node, degree in self.ppi_graph.degree()}

    @staticmethod
    def graph_setup(ppi_file, limit=0.7):
        edge_list_df = pd.read_csv(ppi_file, sep=" ")
        edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
        edge_list_df = edge_list_df[edge_list_df["combined_score"] >= limit]
        G = nx.from_pandas_edgelist(
            edge_list_df,
            "protein1",
            "protein2",
            edge_attr=True
        )
        return G

    def get_degree_dict(self):
        degree_dict = {}
        for node, c_degree in self.ppi_graph.degree():
            if c_degree in degree_dict:
                degree_dict[c_degree].append(node)
            else:
                degree_dict[c_degree] = [node]
        return degree_dict

    def calculate_closest_distances(self, genes_from, genes_to):
        distances = np.full(
            (len(genes_from), len(genes_to)),
            np.inf
        )
        for i, gene_from in enumerate(genes_from):
            for j, gene_to in enumerate(genes_to):
                try:
                    distances[i, j] = nx.shortest_path_length(
                        self.ppi_graph,
                        source=gene_from,
                        target=gene_to
                    )
                except nx.NetworkXNoPath as e:
                    pass
                    # print(e)
        distances = distances.min(axis=1)
        if sum(distances == np.inf) == len(distances):
            dc_sum = np.inf
        else:
            distances = distances[distances != np.inf]
            dc_sum = distances.sum()
        return dc_sum, len(distances)

    def estimate_mu_sigma(self, genes_from, genes_to, degree_diff=.5, samplings=100):
        print("Sampling equivalent sets")

        def _sample_degree_equivalence(gene_set):
            sample_set = []
            for gene in gene_set:
                gene_degree = self.node_degree_dict[gene]
                int_diff = int(gene_degree * degree_diff)
                low = gene_degree - int_diff if int_diff < gene_degree else 1
                max_degree = max(self.degree_dict.keys())
                high = gene_degree + int_diff if gene_degree + int_diff <= max_degree else max_degree

                eq_genes = []
                for idx in range(low, high + 1):
                    if idx in self.degree_dict:
                        eq_genes.append(self.degree_dict[idx])
                    # else:
                    #   print(f"No genes with degree: {idx}")

                eq_genes = list(itertools.chain.from_iterable(eq_genes))
                sample_set += random.sample(eq_genes, 1)

            return sample_set

        distances = np.full(
            samplings,
            np.inf
        )

        for i in range(samplings):
            if i % 20 == 0:
                print(f"Sampling: {i} of {samplings}")
            eq_from = _sample_degree_equivalence(genes_from)
            eq_to = _sample_degree_equivalence(genes_to)
            dc_sum, dc_n = self.calculate_closest_distances(eq_from, eq_to)
            distances[i] = dc_sum/dc_n

        print(f"{sum(distances == np.inf)} of {samplings} samplings ended without path between")

        distances = distances[distances != np.inf]  # remove those without a path between
        return distances.sum(), len(distances), distances.std()

    def calculate_z_scores(self, from_genes, to_genes):
        dab_sum, dab_N = self.calculate_closest_distances(from_genes, to_genes)
        d_sum, d_N, s = self.estimate_mu_sigma(from_genes, to_genes)
        u = d_sum/d_N
        z = (dab_sum/dab_N - u) / s

        return dab_sum, dab_N, u, s, z

    def write_results(self, gene_file_a, gene_file_b, a_name, b_name, output_file):
        with open(gene_file_a, "r") as f:
            genes_a = {l.strip() for l in f.readlines()[1:]}
            genes_a = set(genes_a) & self.ppi_graph.nodes()
        with open(gene_file_b, "r") as f:
            genes_b = {l.strip() for l in f.readlines()[1:]}
            genes_b = set(genes_b) & self.ppi_graph.nodes()

        if os.path.exists(output_file):
            w = open(output_file, "a")
        else:
            w = open(output_file, "w")
            w.write("from\tto\tdsum\tdN\tsab\tu\ts\tz\n")

        dab_sum, dab_N, uab, sab, zab = self.calculate_z_scores(genes_a, genes_b)
        dba_sum, dba_N, uba, sba, zba = self.calculate_z_scores(genes_b, genes_a)

        sab = (dab_sum + dba_sum)/(dab_N + dba_N)

        w.write(f"{a_name}\t{b_name}\t{dab_sum}\t{dab_N}\t{sab}\t{uab}\t{sab}\t{zab}\n")
        w.write(f"{b_name}\t{a_name}\t{dba_sum}\t{dba_N}\t{sba}\t{uba}\t{sba}\t{zba}\n")
