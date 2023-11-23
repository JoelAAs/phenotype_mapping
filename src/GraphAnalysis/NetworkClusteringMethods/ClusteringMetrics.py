import numpy as np
import pandas as pd
import math
import networkx as nx
from scipy.sparse.csgraph import laplacian
from scipy.sparse.linalg import eigsh


def silhouette(edge_file, cluster_file, silhouette_output):
    edges_df = pd.read_csv(edge_file, sep="\t")
    cluster_df = pd.read_csv(cluster_file, sep="\t")
    cluster_dict = {}
    with open(cluster_file, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            node, cluster = line.strip().split("\t")
            if cluster in cluster_dict:
                cluster_dict[cluster].append(node)
            else:
                cluster_dict[cluster] = [node]

    G = nx.from_pandas_edgelist(
        edges_df,
        "query",
        "neighborhood",
        edge_attr=True,
        create_using=nx.DiGraph
    )
    node_list = [n for n in G.nodes]
    adj = nx.adjacency_matrix(G, node_list, weight="probability")
    L = laplacian(adj.todense(), normed=True)  # TODO: version missmatch, needed todense on laptop

    clusters = np.zeros(len(node_list))
    a = np.zeros(len(node_list))
    b = np.zeros(len(node_list))
    s = np.zeros(len(node_list))

    eigenvalues, eigenvectors = eigsh(L, k=len(node_list) - 1)

    for cluster in cluster_dict:
        cluster_idx = [node in cluster_dict[cluster] for node in node_list]
        cluster_dict[cluster] = cluster_idx

    n_eigen = len(cluster_dict) + 2
    if n_eigen > eigenvectors.shape[1]:
        n_eigen = eigenvectors.shape[1]

    for i, node in enumerate(node_list):
        cluster = str(cluster_df[cluster_df["Node"] == node]["Cluster"].values[0])
        clusters[i] = cluster
        node_eigen = eigenvectors[i, 0:n_eigen]
        idx_a = cluster_dict[cluster].copy()
        idx_a[i] = False
        dist_a = (eigenvectors[idx_a, :n_eigen] - node_eigen) ** 2  # works as eigen is np.array
        a[i] = np.sqrt(dist_a.sum(axis=1)).mean()

        idx_b = [not i for i in idx_a]
        idx_b[i] = False
        dist_b = (eigenvectors[idx_b, :n_eigen] - node_eigen) ** 2  # works as eigen is np.array
        b[i] = np.sqrt(dist_b.sum(axis=1)).mean()

        s[i] = (b[i] - a[i]) / max([a[i], b[i]])

    with open(silhouette_output, "w") as w:
        w.write("node\tcluster\ts\ta\tb\n")
        for node, cluster, si, ai, bi in zip(*[node_list, clusters, s, a, b]):
            w.write(f"{node}\t{cluster}\t{si}\t{ai}\t{bi}" + "\n")


