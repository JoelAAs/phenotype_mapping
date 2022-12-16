import pandas as pd
from itertools import combinations
include: "Chembl/Snakefile_chmbl_prediction.smk"


### Drugs
drugs_df = pd.read_csv("data/Drugs.csv", sep="\t")
drugs  = list(drugs_df["ENG"])
drug_combinations = combinations(drugs, 2)
sorted_drug_combinations = [sorted(comb) for comb in drug_combinations]
combination_input = ["work/chembl/distances/{}_{}_all.csv".format(comb[0], comb[1]) for comb in sorted_drug_combinations]

rule all:
    input:
        combination_input