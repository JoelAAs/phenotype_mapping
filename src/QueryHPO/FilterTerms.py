from collections import Counter
import os, glob
import pandas as pd

med_hpo_files = glob.glob("data/micromedex/*txt")

hpo_terms = []
for fn in med_hpo_files:
    with open(fn, "r") as f:
        header = True
        for line in f:
            if header:
                header = False
            else:
                hpo_terms.append(line.strip())

count = Counter(hpo_terms)
df_count = pd.DataFrame(count.items())
df_count = df_count.rename(columns={0: "Term", 1: "N"})
df_terms = pd.read_csv("data/micromedex/all_unique_terms.csv", sep="\t")
df_terms = df_terms.merge(df_count, on=["Term", ], how="right")
df_terms = df_terms[~df_terms.HPO.isnull()]
df_terms = df_terms[df_terms.N > 1]

passed_terms = list(df_terms.HPO)
all_genes = []
empty = []
for term in passed_terms:
    name_term = term.replace(":", "-")
    with open(f"input/HPO/{name_term}.csv", "r") as f:
        genes = [l.strip() for l in f]
        if len(genes) == 1: ## only header
            empty.append(term)
        else:
            all_genes += genes

passed_terms = [t for t in passed_terms if t not in empty]
output_path = "input/HPO-pruned"
gene_count = Counter(all_genes)

for term in passed_terms:
    name_term = term.replace(":", "-")
    with open(f"{output_path}/{name_term}.csv", "w") as w:
        with open(f"input/HPO/{name_term}.csv", "r") as f:
            header = True
            for line in f:
                if header:
                    header = False
                    w.write(line)
                elif gene_count[line.strip()] > 1:
                    w.write(line)
