import pandas as pd
import glob
import os


def get_expected_edgelists(wc):
    edgelists = config["lists"]
    joining_layer = config["joining"]
    expected_lists = [
        "work/edgelists/normalized/{project}_{cutoff}.csv".format(
            project = edgelist,
            cutoff = wc.cutoff
        ) for edgelist in edgelists]
    expected_lists.append(joining_layer)
    return expected_lists

rule filter:
    input:
        edgelist = "work/{project}/term_distance/all_term_distances.csv"
    output:
        filterd = "work/edgelists/filter/{project}_{cutoff}.csv"
    run:
        edgelist = pd.read_csv(input.edgelist, sep = "\t")
        edgelist_filtered = edgelist.query(f"score >= score.quantile(q={wildcards.cutoff})")
        edgelist_filtered.to_csv(output.filterd, sep = "\t", index=False)


rule normalize:
    input:
        filterd = "work/edgelists/filter/{project}_{cutoff}.csv"
    output:
        normalized = "work/edgelists/normalized/{project}_{cutoff}.csv"
    run:
        edgelist_filtered = pd.read_csv(input.filterd, sep="\t")
        edgelist_filtered.score = edgelist_filtered.score/max(edgelist_filtered.score)
        edgelist_filtered.to_csv(output.normalized, sep="\t", index=False)


rule join:
    input:
        get_expected_edgelists
    output:
        "work/edgelists/joined/{joinedname}_{cutoff}.csv"
    shell:
        """
        awk 'FNR > 1' {input} >> {output}
        """
