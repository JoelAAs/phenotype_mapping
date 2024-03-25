from SCnormLapacian import sc_norm_lapacian

checkpoint SCnormLapacian:
    input:
        adj_list = "work/{project}/term_to_term_probability_matrix.csv"
    output:
        clusters = "work/{project}/clustering/SCnorm_{n_clusters}.csv"
    run:
        sc_norm_lapacian(input.adj_list, output.clusters, n_clusters=int(wildcards.n_clusters))

