from ClusteringMetrics import silhouette


rule silhouette_analysis:
    input:
        edge_file="work/{project}/term_to_term_probability_matrix.csv",
        cluster_file = "work/{project}/clustering/SCnorm_{n_clusters}.csv"
    output:
        silhouette_scores = "work/{project}/clustering/metrics/silouette_{n_clusters}.csv"
    run:
        silhouette(input.edge_file, input.cluster_file, output.silhouette_scores)
