"work/{project}/term_gene_distance/full_matrix_{term2}_in_{term1}.csv"

import pandas as pd
import numpy as np
from itertools import permutations

def get_term_permutations(terms):
    return [
        f"work/{{project}}/term_gene_distance/full_matrix_{term2}_in_{term1}.csv"
        for term1, term2 in permutations(terms, 2)
    ]

rule aggregate_Term_probabilities:
    input:
        term_probability_matrixes = get_term_permutations(config["terms"])
    output:
        term_term_probability_matrix = "work/{project}/term_to_term_probability_matrix.csv"
    run:
        terms = config["terms"]
        output_rows = []

        for term_pair_file in input.term_probability_matrixes:
            name_str = term_pair_file.split("_")
            query_set = name_str[-3]
            neighborhood_set = name_str[-1].replace(".csv", "")
            term_pair_df = pd.read_csv(term_pair_file, sep = "\t")
            term_pair_df = term_pair_df.drop('gene', axis=1)
            prob_term2term = term_pair_df.sum(axis=0).mean()
            output_rows.append({
                "query": query_set,
                "neighborhood": neighborhood_set,
                "probability": prob_term2term
            })

        term_prob_df = pd.DataFrame(output_rows)
        term_prob_df.to_csv(output.term_term_probability_matrix, sep = "\t")