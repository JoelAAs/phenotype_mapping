import pandas as pd
from PermutationAnalysis import sc_norm_permutate, calculate_correct_node_ratio, score_hpo_terms


# TODO: Make hpo node scoring based on connected nodes

rule permutation_results:
    params:
        n_clusters = 3,
        n_permutations = 1000,
        fraction_forgotten = 0.1
    input:
        edge_file = "work/{project}/term_to_term_probability_matrix.csv"
    output:
        permutations = "work/{project}/clustering/metrics/permutations.csv"
    run:
        permutation_results_df = sc_norm_permutate(
            input.edge_file,
            n_clusters=params.n_clusters,
            n_permut=params.n_permutations,
            fraction_forgotten=params.fraction_forgotten
        )
        permutation_results_df.to_csv(output.permutations, sep="\t", index=False)

rule group_permutation_robustness:
    input:
        edge_file = "work/{project}/term_to_term_probability_matrix.csv",
        node_positions = "data/Node_groups.csv",
        permutations= "work/{project}/clustering/metrics/permutations.csv"
    output:
        clustering_ratio = "work/{project}/clustering/metrics/sensitivity.csv"
    run:
        permutation_results_df = pd.read_csv(input.permutations, sep = "\t")
        node_vote_group_df = calculate_correct_node_ratio(
            permutation_results=permutation_results_df,
            node_grouping_file=input.node_positions
        )
        sensitivity_series = node_vote_group_df.mean(skipna=True, axis=0)

        with open(output.clustering_ratio, "w") as w:
            w.write("Node\tSensitivity\n")
            for node, sens in sensitivity_series.items():
                w.write(f"{node}\t{sens}\n")


rule score_HPO:
    params:
        n_clusters = 3,
        n_permutations = 1000,
        fraction_forgotten = 0.1
    input:
        permutations = "work/{project}/clustering/metrics/permutations.csv",
        node_positions = "data/Node_groups.csv",
        drug_adr_csv = "data/hpo_pruned_med_directed.csv"
    output:
        senes_hpo = "work/{project}/clustering/metrics/hpo_sensitivity.csv"
    run:
        hpo_sensitivity_series = score_hpo_terms(input.permutations, input.node_positions, input.drug_adr_csv)
        with open(output.senes_hpo, "w") as w:
            w.write("Node\tSensitivity\n")
            for node, sens in hpo_sensitivity_series.items():
                w.write(f"{node}\t{sens}\n")