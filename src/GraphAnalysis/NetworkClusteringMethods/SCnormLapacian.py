import pandas as pd
import networkx as nx
from matplotlib.pylab import cm
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.cluster import KMeans
from scipy.sparse.csgraph import laplacian
from scipy.sparse.linalg import eigsh


def sc_norm_lapacian(edge_file, output, n_clusters=4):
    edge_list_df = pd.read_csv(edge_file, sep="\t")
    G = nx.from_pandas_edgelist(
        edge_list_df,
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

    partition_dict = {k: [] for k in range(n_clusters)}
    for i, node in enumerate(node_list):
        partition_dict[km.labels_[i]].append(node)

    with open(output, "w") as w:
        w.write("Node\tCluster\n")
        for cluster in partition_dict:
            for node in partition_dict[cluster]:
                w.write(f"{node}\t{cluster}\n")
