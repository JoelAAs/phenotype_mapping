import pandas as pd
import bz2
import glob
import os
from itertools import combinations

### CONFIG
input_files = glob.glob(f"{config['input_location']}/{config['project']}/*csv")
term_combinations = [
    os.path.basename(
        term
    ).replace(".csv", "") for term in input_files
]
term_combinations = combinations(term_combinations, 2)


### RULE
rule get_term_interaction_gene_interaction:
    input:
        term_a = "{input_location}/{{project}}/{{term_a}}.csv".format(input_location=config["input_location"]),
        term_b= "{input_location}/{{project}}/{{term_b}}.csv".format(input_location=config["input_location"])
    output:
        gene_interactions = "work/{project}/gene_interactions/{term_a}_{term_b}.csv.bz2"
    run:
        df_a = pd.read_csv(input.term_a, sep="\t")
        df_b = pd.read_csv(input.term_b, sep="\t")
        df_combinations = df_a.merge(df_b, how="cross", suffixes=("_a", "_b"))
        df_combinations.to_csv(output.gene_interactions, sep = "\t", index=False)


rule get_unique_gene_combinations:
    input:
        all_combination_dfs = [f"work/{{project}}/gene_interactions/{term_a}_{term_b}.csv.bz2"
                   for term_a, term_b in term_combinations]
    output:
        unique_gene_pairs = "work/{project}/gene_interactions/unique.csv",
        all_unique_genes = "work/{project}/gene_interactions/unique_genes.csv"
    run:
        unique_set = set()
        unique_genes = set()
        for interaction_file_name in input.all_combination_dfs:
            with bz2.open(interaction_file_name,"r") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.decode("utf-8").strip().split("\t")
                    gene_a, gene_b = line[:2]
                    pair = tuple([line[0], line[1]])
                    pair_rev = tuple([line[1], line[0]])
                    if gene_a not in unique_genes:
                        unique_genes.update({gene_a})
                    if gene_b not in unique_genes:
                        unique_genes.update({gene_b})
                    if pair not in unique_set:
                        unique_set.update({pair})
                    if pair_rev not in unique_set:
                        unique_set.update({pair_rev})

        with open(output.unique_gene_pairs, "w") as w:
            w.write("gene_a\tgene_b\n")
            for a, b in unique_set:
                w.write(f"{a}\t{b}\n")

        with open(output.all_unique_genes, "w") as w:
            w.write("gene\n")
            for gene in unique_genes:
                w.write(f"{gene}\n")
