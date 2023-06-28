import pandas as pd
import markov_clustering as mc
import networkx as nx
from matplotlib.pylab import show, cm, axis

edge_list_df = pd.read_csv("work/edgelists/joined/test_0.9.csv", sep="\t")
# edge_list_df = edge_list_df[edge_list_df.score > edge_list_df.score.mean()]
G = nx.from_pandas_edgelist(edge_list_df, "term1", "term2", edge_attr=True)

labels = list(G.nodes())
matrix = nx.adjacency_matrix(G, weight="score")
result = mc.run_mcl(matrix, inflation=1.3)
clusters = mc.get_clusters(result)

# map node to cluster id for colors
cluster_map = {node: i for i, cluster in enumerate(clusters) for node in cluster}
colors = [cluster_map[i] for i in range(len(G.nodes()))]

nx.draw_networkx(G,
                 node_color=colors,
                 node_size=50,
                 with_labels=True,
                 edge_color="silver",
                 cmap=cm.tab20)

axis("off")
show(block=False)
