import pandas as pd
import networkx as nx
import numpy as np
from sklearn.cluster import KMeans
from scipy.sparse.csgraph import laplacian
from scipy.sparse.linalg import eigsh
import random
from collections import Counter


def sc_norm_permutate(edge_file, n_clusters=3, n_permut=1000, fraction_forgotten=.1):
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
        if i % 50 == 0:
            print(f"{i}/{n_permut} done")
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
    group_dict = dict()
    size_dict = dict()
    with open(node_grouping_file, "r") as f:
        lines = [l.strip() for l in f]
        for line in lines[1:]:
            node, group = line.split("\t")
            group_dict[node] = group
            if group not in size_dict:
                size_dict[group] = 1
            else:
                size_dict[group] += 1

    column_to_keep = list(group_dict.keys())
    permutation_results = permutation_results[column_to_keep]
    correct_votes = pd.DataFrame(
        np.full(permutation_results.shape, np.nan),
        columns=permutation_results.columns.values
    )

    for i, row in permutation_results.iterrows():
        row_cluster = dict()
        row_voting = dict()
        for node, cluster in row.items():
            if cluster not in row_voting:
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
                correct_votes.at[i, node] = int(row_cluster[cluster] == group_dict[node])


    return correct_votes


def score_hpo_terms(permutation_csv, node_positions, drug_adr_csv):
    permutation_df = pd.read_csv(permutation_csv, sep="\t").T
    node_positions_df = pd.read_csv(node_positions, sep="\t")
    drug_adr_df = pd.read_csv(drug_adr_csv, sep="\t")
    drug_adr_df = drug_adr_df[~drug_adr_df["from"].isin(node_positions_df["Node"])]
    adrs = set(drug_adr_df["from"])

    percent_hits = pd.DataFrame(
        np.full((permutation_df.shape[1], len(adrs)), np.nan),
        columns=adrs
    )

    for j in range(permutation_df.shape[1]):
        for adr_node in adrs:
            med_connect_adr = drug_adr_df[drug_adr_df["from"] == adr_node]["to"]
            chosen_cluster = permutation_df.loc[adr_node][j]
            terms_in_cluster = permutation_df[permutation_df[j] == chosen_cluster][j].index.tolist()
            percent_associated = sum([True for m in med_connect_adr if m in terms_in_cluster])/len(med_connect_adr)
            percent_hits.at[j, adr_node] = percent_associated

    return percent_hits.mean()