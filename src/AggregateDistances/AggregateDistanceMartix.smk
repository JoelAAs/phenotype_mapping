import os
import glob
from itertools import combinations
import pandas as pd

input_files = glob.glob(f"{config['input_location']}/{config['project']}/*csv")
term_combinations = [
    os.path.basename(
        term
    ).replace(".csv", "") for term in input_files
]
term_combinations = list(combinations(term_combinations, 2))

rule get_term_distance_matrix:
    input:
        all_combinations = [
            f"work/{{project}}/gene_distance/{term_a}_{term_b}.csv.bz2" for
            term_a, term_b in term_combinations]
    output:
        all_distance = "work/{project}/term_distance/all_term_distances.csv"
    run:
        with open(output.all_distance, "w") as w:

            for term_a, term_b in term_combinations:
                df_distance = pd.read_csv(
                    f"work/{wildcards.project}/gene_distance/{term_a}_{term_b}.csv.bz2", sep = "\t"
                )
                print(df_distance)
                tot_sum = sum(df_distance["weight_a"].astype(float) + df_distance["weight_b"].astype(float))
                n_connections = len(df_distance)
                avg_weight = tot_sum/n_connections
                w.write(f"{term_a}\t{term_b}\t{tot_sum}\t{tot_sum}\t{n_connections}\t{avg_weight}\n")
