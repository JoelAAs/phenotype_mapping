import pandas as pd

from SeparationDistance import CalculateSeparation
import glob
from itertools import combinations
import os
from scipy.stats import norm

n_batches = 8

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


rule batch_combinations:
    input:
        terms= lambda wc: get_benchmark_input(wc)
    output:
        batches = expand("work/{{project}}/benchmarking/batch/combinations_{n}.csv", n=range(n_batches))
    run:
        combs = list(combinations(input.terms,2))
        steps = int(len(combs)/n_batches)
        i = 0
        b = 0
        w = open(output.batches[b], "w")
        for comb in combs:
            if i > steps:
                w.close()
                b += 1
                w = open(output.batches[b], "w")
                i = 0
            w.write("\t".join(comb) + "\n")
            i +=1
        w.close()




rule generate_z_scores:
    params:
        ppi_file = "data/9606.protein.links_above_700.v11.5.txt"
    input:
        input_batch = "work/{project}/benchmarking/batch/combinations_{n}.csv"
    output:
        z_scores = "work/{project}/benchmarking/z_score/Affinity_zscore_batch_{n}.csv"
    run:
        with open(input.input_batch, "r") as f:
            for line in f:
                comb = line.strip().split("\t")

                name_a = os.path.basename(comb[0]).replace(".csv", "")
                name_b = os.path.basename(comb[1]).replace(".csv","")
                cs = CalculateSeparation(params.ppi_file)
                cs.write_results(comb[0], comb[1], name_a, name_b, output.z_scores)

rule aggregate_results:
    input:
        result_batches = expand("work/{{project}}/benchmarking/z_score/Affinity_zscore_batch_{n}.csv", n=range(n_batches))
    output:
        z_scores = "work/{project}/benchmarking/z_score/Affinity_zscore.csv"
    shell:
        """
        echo "from\tto\tdsum\tdN\tsab\tu\ts\tz\n" > {output.z_scores}
        cat {input.result_batches} >> {output.z_scores}
        """


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
        max_sab = sab_df["sab"].max()
        sab_df["sab"] = sab_df["sab"].apply(lambda x: 1 - x/max_sab) # Transform to affinity for spectral clustering
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