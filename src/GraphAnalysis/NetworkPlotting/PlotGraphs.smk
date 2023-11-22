from plot_graph import graph_plot, cluster_graph_plot
import pandas as pd
import os
import glob

def get_clusteredge_output(wc):
    cluster_edge_ck = checkpoints.ClusterEdges.get(**wc).output[0]
    clusters, = glob_wildcards(os.path.join(cluster_edge_ck, "SC_edges_{cluster}.csv"))
    return expand(os.path.join(cluster_edge_ck, "SC_edges_{cluster}.csv"), cluster=clusters)

rule SCPlot:
    input:
        graph = "work/{project}/term_to_term_probability_matrix.csv",
        clusters = "work/{project}/clustering/SCnorm.csv"
    output:
        figure = "work/{project}/clustering/plots/SC.png"
    run:
        graph_plot(
            input.graph,
            input.clusters,
            wildcards.project,
            "Spectral Clustering",
            output.figure)

checkpoint ClusterEdges:
    input:
        graph="work/{project}/term_to_term_probability_matrix.csv",
        clusters="work/{project}/clustering/SCnorm.csv"
    output:
        directory("work/{project}/clustering/clustered_edges")
    run:
        os.mkdir(output[0])

        cluster_dict = {}
        with open(input.clusters) as f:
            lines = f.readlines()
            for line in lines[1:]:
                node, cluster = line.strip().split("\t")
                if cluster in cluster_dict:
                    cluster_dict[cluster].append(node)
                else:
                    cluster_dict[cluster] = [node]

        edge_df = pd.read_csv(input.graph, sep = "\t")
        for cluster in cluster_dict:
            current_nodes = cluster_dict[cluster]

            cluster_df = edge_df[edge_df["query"].isin(current_nodes)]
            cluster_df = cluster_df[cluster_df["neighborhood"].isin(current_nodes)]
            cluster_df.to_csv(
                f"work/{wildcards.project}/clustering/clustered_edges/SC_edges_{cluster}.csv",
                sep="\t",
                index=None,
            )

rule ClusterPlot:
    input:
        get_clusteredge_output
    output:
        check = "work/{project}/clustering/plots/single_clusters/done.txt"
    run:
        with open(output.check, "w") as w:
            w.write("check")
        for i, cluster_edge_file in enumerate(input):
            print(input)
            cluster_graph_plot(
                cluster_edge_file,
                wildcards.project,
                "Spectral clustering",
                f"work/{wildcards.project}/clustering/plots/single_clusters/cluster_{i}.png"
            )
