import pandas as pd
import random
import numpy as np


class GraphKNN:
    def __init__(self, adjacency_list_filename, k=3, iterations=10):
        self.iterations = iterations
        self.adj_df = pd.read_csv(adjacency_list_filename, sep="\t")
        nodes = list(set(self.adj_df["query"]))
        n_per_cluster = round(len(nodes)/k)
        clusters = np.array([0]*len(nodes))
        for i in range(k):
            if i == k-1:
                clusters[i*n_per_cluster:] = i
            else:
                clusters[i*n_per_cluster:(i+1)*n_per_cluster] = i

        random.shuffle(clusters)

        d = {
            "neighborhood": nodes,
            "cluster": clusters
        }
        self.cluster_df = pd.DataFrame(d)

    def step(self):
        cluster_edges = self.adj_df.merge(
            self.cluster_df, on="neighborhood"
        )
        sum_prob = cluster_edges.groupby(["query", "cluster"], as_index=False).sum("probability")
        sum_prob["max_val"] = sum_prob.groupby(["query"], as_index=False)["probability"].transform(max)
        sum_prob = sum_prob[sum_prob["max_val"] == sum_prob["probability"]][["query", "cluster"]]
        self.cluster_df = sum_prob.rename({
            "query": "neighborhood",
            "cluster": "cluster"
        }, axis=1)

    def cluster(self, output_cluster_location):
        for _ in range(self.iterations):
            self.step()

        out_df = self.cluster_df.rename({
            "neighborhood": "node",
            "cluster": "cluster"
        }, axis=1)
        out_df.to_csv(output_cluster_location, sep="\t", index_label=False)
