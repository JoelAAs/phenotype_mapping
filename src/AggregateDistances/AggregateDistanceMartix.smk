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
            f"work/{{project}}/gene_distance/{term_a}_{term_b}.csv" for
            term_a, term_b in term_combinations]
    output:
        all_distance = "work/{project}/term_distance/all_term_distances.csv"
    run:
        with open(output.all_distance, "w") as w:
            w.write("term_a\tterm_2\tscore\n")
            for term_a, term_b in term_combinations:
                df_distance = pd.read_csv(
                    f"work/{wildcards.project}/gene_distance/{term_a}_{term_b}.csv", sep = "\t"
                )
                score = df_distance.groupby("from_term")["probability"].sum().sum()/2
                w.write(f"{term_a}\t{term_b}\t{score}\n")
