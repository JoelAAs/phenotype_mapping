from PermutationAnalysis import sc_norm_permutate, calculate_correct_node_ratio


# TODO: Sort out all that exists within cluster
# TODO: Make hpo node scoring based on connected nodes
rule node_permutation_robustness:
    params:
        n_clusters = 3,
        n_permutations = 1000,
        fraction_forgotten = 0.1
    input:
        edge_file = "work/{project}/term_to_term_probability_matrix.csv",
        node_positions = "data/Node_groups.csv"
    output:
        clustering_ratio = "work/{project}/clustering/metrics/sensitivity.csv"
    run:
        permutation_results_df = sc_norm_permutate(
            input.edge_file,
            n_clusters=params.n_clusters,
            n_permut=params.n_permutations,
            fraction_forgotten=params.fraction_forgotten
        )
        node_vote_group_df = calculate_correct_node_ratio(
            permutation_results=permutation_results_df,
            node_grouping_file=input.node_positions
        )
        sensitivity_series = node_vote_group_df.mean(skipna=True, axis=1)

        with open(output.clustering_ratio, "w") as w:
            w.write("Node\tSensitivity\n")
            for node, sens in sensitivity_series.items():
                w.write(f"{node}\t{sens}\n")