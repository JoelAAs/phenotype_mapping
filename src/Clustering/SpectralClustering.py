from sklearn.cluster import SpectralClustering
import pandas as pd
import networkx as nx
from matplotlib.pylab import show, cm, axis


edge_list_df = pd.read_csv("work/edgelists/joined/test_0.9.csv", sep="\t")

G = nx.from_pandas_edgelist(edge_list_df, "term1", "term2", edge_attr=True)
adj_matrix = nx.adjacency_matrix(G, weight=None)

spec_cl = SpectralClustering(
    3,
    affinity='precomputed',
    n_init=100
)
spec_cl.fit(adj_matrix)

nx.draw_networkx(
    G,
     node_color=spec_cl.labels_,
     node_size=50,
     with_labels=True,
     edge_color="silver",
     cmap=cm.tab20)

axis("off")
#plt.savefig("work/edgelists/plots/plot.png")
show(block=False)
