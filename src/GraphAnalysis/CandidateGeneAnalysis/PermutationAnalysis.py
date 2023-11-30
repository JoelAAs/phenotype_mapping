import pandas as pd
import networkx as nx
import numpy as np
from sklearn.cluster import KMeans
from scipy.sparse.csgraph import laplacian
from scipy.sparse.linalg import eigsh
import random
from collections import Counter


def sc_norm_permutate(edge_file, n_clusters=4, n_permut=1000, fraction_forgotten=.1):
    edge_list_df = pd.read_csv(edge_file, sep="\t")
    nodes = list(set(
        edge_list_df["query"].tolist() +
        edge_list_df["neighborhood"].tolist()
    ))
    clusters_mat = np.full((n_permut, len(nodes)), np.nan)
    permutation_results_df = pd.DataFrame(
        clusters_mat,
        columns=nodes
    )

    for i in range(n_permut):
        sample_df = edge_list_df.sample(frac=1 - fraction_forgotten)

        G = nx.from_pandas_edgelist(
            sample_df,
            "query",
            "neighborhood",
            edge_attr=True,
            create_using=nx.DiGraph
        )
        node_list = [n for n in G.nodes]
        adj = nx.adjacency_matrix(G, node_list, weight="probability")
        L = laplacian(adj.todense(), normed=True, use_out_degree=True)  # in as default

        eigenvalues, eigenvectors = eigsh(L, k=len(node_list) - 1)
        km = KMeans(
            n_clusters=n_clusters,
            n_init=60,
            max_iter=1000
        ).fit(
            eigenvectors[:, :n_clusters]
        )

        for j, node in enumerate(node_list):
            permutation_results_df.at[i, node] = km.labels_[j]

    return permutation_results_df


def calculate_correct_node_ratio(permutation_results, node_grouping_file):
    group_dict = {}
    size_dict = {}
    with open(node_grouping_file, "r") as f:
        lines = [l.strip() for l in f]
        for line in lines[1:]:
            node, group = line.split("\t")
            group_dict[node] = group
            if group in size_dict:
                size_dict[group] = 1
            else:
                size_dict[group] += 1

    correct_votes = pd.DataFrame(
        np.full(permutation_results.shape, np.nan),
        columns=permutation_results.columns.values
    )

    row_cluster = {}
    row_voting = {}
    for i, row in permutation_results.iterrows():
        for node, cluster in row.items():
            if cluster in row_voting:
                row_voting[cluster] = [group_dict[node]]
            else:
                row_voting[cluster].append(group_dict[node])

        for cluster in row_voting:
            votes = Counter(row_voting[cluster])
            current_max = 0
            max_group = ""
            for group in votes:
                size_weighted_votes = votes[group] / size_dict[group]
                if size_weighted_votes >= current_max:
                    current_max = size_weighted_votes
                    max_group = group
            row_cluster[cluster] = max_group

        for node, cluster in row.items():
            if cluster is not np.nan:
                correct_votes.at[i, node] = row_cluster[cluster] == group_dict[node]

    return correct_votes
