import pandas as pd
import networkx as nx
from matplotlib.pylab import cm
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

sns_cm = sns.color_palette("tab10").as_hex()
width = 7.3
height = 7.3
fig, axes = plt.subplots(1, 1, figsize=(width, height), layout="tight")

cutoff = 0.80

edges = f"work/edgelists/joined/full-drugbank_{cutoff}.csv"
clusters = f"work/edgelists/clustering/table/full-drugbank_{cutoff}_cutoffMCL.csv"
edge_list_df = pd.read_csv(edges, sep="\t", header=None)
edge_list_df = edge_list_df.rename({0: "term1", 1: "term2", 2: "score"}, axis=1)
G = nx.from_pandas_edgelist(edge_list_df, "term1", "term2", edge_attr=True)

label_df = pd.read_csv(clusters, sep="\t")
lab_dict = dict()
for i, row in label_df.iterrows():
    lab_dict[row[0]] = row[1]

for node in G:
    G.nodes[node]["shape"] = ('o' if node[0:2] in ["HP", "OR"] else "d")
    G.nodes[node]["size"] = (100 if node[0:2] == "HP" else 150)
    G.nodes[node]["color"] = sns_cm[lab_dict.get(node)]

# plt.savefig(f"paper-prep/plots/MLC_{cutoff}.png", dpi=500)
pos = nx.spring_layout(G)
node_shapes = ["o", "d"]

nx.draw_networkx_edges(G, pos, edge_color="grey")  # draw edges

# Draw the nodes for each shape with the shape specified
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

plt.title("MCL, 1.6 inflation at 20% highest weighted edges")
plt.savefig("paper-prep/plots/MLC0.2.png", dpi=300)
