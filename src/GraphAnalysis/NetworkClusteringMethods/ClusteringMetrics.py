import numpy as np
import pandas as pd
import math

def silhouette(edge_file, cluster_file):
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

    edges_df = edges_df.merge(cluster_file, right_on="query", left_on="Node")
    edges_df["distance"] = 1 - edges_df["probability"]

    nodes = list(set(edges_df["query"]))
    clusters = np.zeros(len(nodes))
    a = np.zeros(len(nodes))
    b = np.zeros(len(nodes))
    s = np.zeros(len(nodes))

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
    L = laplacian(adj, normed=True)  # in as default

    for i, node in enumerate(nodes):
        cluster = cluster_df[cluster_df["Node"] == node]["Cluster"][0]
        clusters[i] = cluster
        a_list = edges_df[edges_df["Cluster"] == cluster][edges_df["query"] == node]["distance"].tolist()
        a_list += [1]*(len(cluster_dict[cluster]) - 1)
        a[i] = sum(a_list)/len(a_list)
        b_list = edges_df[edges_df["Cluster"] != cluster][edges_df["query"] == node]["distance"].tolist()
        b_list += [1] * (len(nodes) - len(cluster_dict[cluster]))
        a[i] = sum(b_list) / len(b_list)
        s[i] = (b[i] - a[i]) / max([a[i], b[i]])

    for node, cluster, si, ai, bi in zip([nodes, clusters, s, a, b]):

