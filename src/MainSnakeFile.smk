import pandas as pd
from itertools import combinations
include: "Chembl/Snakefile_chmbl_prediction.smk"


### Drugs
def  get_all_medication(file_csv):
    df = pd.read_csv(file_csv, sep = "\t")
    return list(df.ENG)

drugs = get_all_medication("data/Drugs.csv")
print(f"there are {len(drugs)} drugs.")

drug_combinations = combinations(drugs, 2)
sorted_drug_combinations = [sorted(comb) for comb in drug_combinations]
combination_input = ["work/chembl/distances/{}_{}_all.csv.bz2".format(comb[0], comb[1]) for comb in sorted_drug_combinations]


print(f"and {len(combination_input)} combinations.")

rule all:
    input:
        #"work/chembl/pairs/distance_pairs_possible_ends.csv"
        #expand("work/chembl/predicted_interaction/annotated/{drug}.csv", drug = drugs),
        #combination_input
        "work/chembl/distances/summed_scores.csv"
