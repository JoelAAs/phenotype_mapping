import pandas as pd

from SeparationDistance import CalculateSeparation
import glob
from itertools import combinations
import os
from scipy.stats import norm

projects = [
    #"HPO-pruned",
    "full-drugbank-benchmark",
    #"MedAdr"
]

def get_benchmark_input(wc):
    input_folder = {
        "full-drugbank-benchmark": "full-drugbank",
        "HPO-pruned-benchmark": "HPO-pruned"
    }
    return glob.glob(f"input/{input_folder[wc.project]}/*.csv")


rule generate_z_scores:
    params:
        ppi_file = "data/9606.protein.links_above_700.v11.5.txt"
    input:
        terms = lambda wc: get_benchmark_input(wc)
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


rule translate_as_pd:
    input:
        z_scores = "work/{project}/benchmarking/z_score/Affinity_zscore.csv"
    output:
        edges_z = "work/{project}/term_to_term_probability_matrix_zscore.csv",
        edges_sab = "work/{project}/term_to_term_probability_matrix_sab.csv"
    run:
        z_df = pd.read_csv(input.z_scores, sep="\t")
        z_df["probability"] = z_df.z.apply(lambda x: 1 - norm.cdf(x,0,1))
        sab_df = z_df[["from", "to", "sab"]]
        mean_sab = sab_df["sab"].mean()
        sab_df["sab"] = sab_df["sab"].apply(lambda x: 2*mean_sab - x) # flip distribution => high value mean similar
        sab_df = sab_df[["from", "to", "sab"]].rename(
            {
                "from": "query",
                "to": "neighborhood",
                "sab": "sab"}, axis=1
        )

        sab_df.to_csv(output.edges_sab, sep = "\t", index=None)

        z_df = z_df[["from", "to", "probability"]].rename(
            {
                "from" :"query",
                "to":"neighborhood",
                "probability": "probability"}, axis=1
        )
        z_df.to_csv(output.edges_z, sep = "\t", index=None)