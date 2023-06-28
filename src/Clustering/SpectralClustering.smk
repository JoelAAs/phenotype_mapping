from sklearn.cluster import SpectralClustering
import pandas as pd
import networkx as nx
from matplotlib.pylab import cm
import matplotlib.pyplot as plt



rule SpectralClustering:
    params:
        n_cluster = config["number_clusters"]
    input:
        edges = "work/edgelists/joined/{joinedname}_{cutoff}.csv"
    output:
        plot = "work/edgelists/clustering/plots/{joinedname}_{cutoff}_SC.png",
        table = "work/edgelists/clustering/table/{joinedname}_{cutoff}_SC.csv"
    run:
        edge_list_df = pd.read_csv(input.edges, sep="\t", header=None)
        edge_list_df = edge_list_df.rename({0: "term1", 1: "term2", 2: "score"}, axis = 1)
        G = nx.from_pandas_edgelist(edge_list_df, "term1", "term2", edge_attr=True)
        adj_matrix = nx.adjacency_matrix(G, weight=None)
        spec_cl = SpectralClustering(
            params.n_cluster,
            affinity='precomputed',
            n_init=100
        )
        spec_cl.fit(adj_matrix)
        nx.draw_networkx(
            G,
            node_color=spec_cl.labels_,
            font_size = 5,
            node_size=50,
            with_labels=True,
            edge_color="silver",
            cmap=cm.tab20)
        plt.savefig(output.plot, dpi=500)

        with open(output.table, "w") as w:
            w.write("term\tcluster\n")
            for label, cluster in zip(list(G.nodes), spec_cl.labels_):
                w.write(f"{label}\t{cluster}\n")


