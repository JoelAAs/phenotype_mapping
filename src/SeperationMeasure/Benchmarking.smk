from SeperationDistance import CalculateSeparation
import glob
from itertools import combinations
import os

projects = [
    "HPO-pruned",
    "full-drugbank",
    #"MedAdr"
]

def get_term_input(wc):
    return glob.glob(f"input/{wc.project}/*.csv")

rule all:
    input:
        expand("work/{project}/benchmarking/z_score/Affinity_zscore.csv", project = projects)


rule generate_z_scores:
    params:
        ppi_file = "data/9606.protein.links_above_700.v11.5.txt"
    input:
        terms = lambda wc: get_term_input(wc)
    output:
        z_scores = "work/{project}/benchmarking/z_score/Affinity_zscore.csv"
    run:
        print(input.terms)
        i = 0
        combs = list(combinations(input.terms,2))
        for comb in combs:
            name_a = os.path.basename(comb[0]).replace(".csv", "")
            name_b = os.path.basename(comb[1]).replace(".csv","")
            print(f"checking combination {i} of {len(combs)}")
            i += 1
            cs = CalculateSeparation(params.ppi_file)
            cs.write_results(comb[0], comb[1], name_a, name_b, output.z_scores)