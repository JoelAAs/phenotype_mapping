from GraphKNN import GraphKNN
from SCnormLapacian import sc_norm_lapacian

rule SCnormLapacian:
    input:
        adj_list = "work/{project}/term_to_term_probability_matrix.csv"
    output:
        clusters = "work/{project}/clustering/SCnorm_{n_clusters}.csv"
    run:
        sc_norm_lapacian(input.adj_list, output.clusters, n_clusters=int(wildcards.n_clusters))

rule clusterKNN:
    input:
        adj_list = "work/{project}/term_to_term_probability_matrix.csv"
    output:
        clusters = "work/{project}/clustering/KNN/Clusters_{iterations}_KNN.csv"
    run:
        gnn = GraphKNN(input.adj_list, iterations=int(wildcards.iterations))
        gnn.cluster(output.clusters)