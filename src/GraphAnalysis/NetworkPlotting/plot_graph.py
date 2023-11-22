import pandas as pd
import networkx as nx
from matplotlib.pylab import cm
import matplotlib.pyplot as plt
import seaborn as sns


def graph_plot(edge_file, cluster_file, project, method, output):
    with open(cluster_file, "r") as f:
        clusters_lines = [line.strip() for line in f.readlines()]

    cluster_dict= {}
    for line in clusters_lines[1:]:
        node, cluster = line.split("\t")
        cluster = int(cluster)
        if cluster in cluster_dict:
            cluster_dict[cluster].append(node)
        else:
            cluster_dict[cluster] = [node]

    edge_list_df = pd.read_csv(edge_file, sep="\t")
    G = nx.from_pandas_edgelist(
        edge_list_df,
        "query",
        "neighborhood",
        edge_attr=True,
        create_using=nx.DiGraph
    )

    sns_cm = sns.color_palette("tab10").as_hex()
    sns.set_theme()
    cluster_colors = sns_cm[:len(cluster_dict)]
    width = 7.3
    height = 7.3
    fig, axes = plt.subplots(1, 1, figsize=(width, height))

    pos = nx.spring_layout(G, weight="probability")
    for cluster in cluster_dict:
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=cluster_dict[cluster],
            node_color=cluster_colors[cluster]
        )
    nx.draw_networkx_edges(G, pos, edge_color="grey")
    labels = {node: node for node in G}
    nx.draw_networkx_labels(G, pos, labels, font_color="black")
    plt.title(f"{method} results of {project}")

    plt.savefig(output, dpi=300)


def cluster_graph_plot(edge_file, project, method,  output):

    edge_list_df = pd.read_csv(edge_file, sep="\t")
    G = nx.from_pandas_edgelist(
        edge_list_df,
        "query",
        "neighborhood",
        edge_attr=True,
        create_using=nx.DiGraph
    )
    sns_cm = sns.color_palette("tab10").as_hex()

    for node in G:
        G.nodes[node]["shape"] = ('o' if node[0:2] in ["HP", "OR"] else "d")
        G.nodes[node]["size"] = (100 if node[0:2] in ["HP", "OR"] else 150)
        G.nodes[node]["color"] = (sns_cm[0] if node[0:2] in ["HP", "OR"] else sns_cm[1])


    width = 7.3
    height = 7.3
    fig, axes = plt.subplots(1, 1, figsize=(width, height))
    node_shapes = ["o", "d"]
    pos = nx.spring_layout(G, weight="probability")

    for shape in set(node_shapes):
        # the nodes with the desired shapes
        node_list = [node for node in G.nodes() if G.nodes[node]['shape'] == shape]
        if shape == "d":
            labels = {node: node for node in G if G.nodes[node]["shape"] == "d"}
            nx.draw_networkx_labels(G, pos, labels, font_color="black")

        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=node_list,
            node_size=[G.nodes[node]['size'] for node in node_list],
            node_color=[G.nodes[node]['color'] for node in node_list],
            node_shape=shape,
            edgecolors="grey"
        )

    nx.draw_networkx_edges(G, pos, edge_color="grey")
    labels = {node: node for node in G}
    plt.title(f"{method} results of {project}, single cluster")

    plt.savefig(output, dpi=300)