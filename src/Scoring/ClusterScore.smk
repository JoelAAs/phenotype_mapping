import pandas as pd
import numpy as np

rule calculate_score:
    input:
        term_clustering = "work/edgelists/clustering/table/{joinedname}_{param}_{method}.csv",
        cluster_identity = "data/scoring/med-cluster.csv",
        hpo_score = "data/scoring/HPO-count-cluster.csv"
    output:
        part = "work/edgelists/clustering/results_parts/{joinedname}_{param}_{method}.csv"
    run:
        cluster_id = pd.read_csv(input.cluster_identity, sep = "\t")
        hpo_score = pd.read_csv(input.hpo_score, sep = "\t")
        hpo_max_score = hpo_score.groupby("HPO").percent_count.max().sum()
        term_clustering = pd.read_csv(input.term_clustering, sep = "\t")

        cluster_id_assigned = cluster_id.merge(
            term_clustering, left_on="drug", right_on="term", how="left", suffixes=["_id", "_assigned"])

        cluster_id_assigned_grouped = cluster_id_assigned.groupby(["cluster_id", "cluster_assigned"], as_index=False).term.count()
        cluster_id_assigned_grouped = cluster_id_assigned_grouped.loc[
            cluster_id_assigned_grouped.groupby("cluster_id").term.idxmax()
        ]
        meds_in_groups = cluster_id.groupby("cluster", as_index=False).count()
        meds_in_groups = meds_in_groups.rename({"cluster": "cluster_id", "drug": "max_n_terms"}, axis = 1)
        cluster_id_assigned_grouped = cluster_id_assigned_grouped.merge(meds_in_groups,
            on="cluster_id", how="right"

        )
        # fist cluster to claim get it
        cluster_id_assigned_grouped.loc\
            [cluster_id_assigned_grouped.cluster_assigned.duplicated(), ["cluster_assigned", "term"]] = None

        cluster_id_assigned_grouped.term = cluster_id_assigned_grouped.term.fillna(0)
        medicine_hitrate = cluster_id_assigned_grouped.term.sum()/cluster_id_assigned_grouped.max_n_terms.sum()

        cluster_id = \
            cluster_id_assigned_grouped[~cluster_id_assigned_grouped.cluster_assigned.isna()][["cluster_id", "cluster_assigned"]]
        cluster_id = cluster_id.rename({"cluster_id": "cluster", "cluster_assigned": "cluster_assigned"}, axis=1)
        hpo_score = hpo_score.merge(cluster_id, on="cluster", how = "left")
        term_clustering = term_clustering.rename({"term": "HPO", "cluster": "cluster_assigned"}, axis = 1)
        hpo_score = hpo_score.merge(term_clustering, on=["HPO", "cluster_assigned"], how="right")
        hpo_summed_score = hpo_score[~hpo_score["count"].isna()]["percent_count"].sum()

        hpo_percent_score = hpo_summed_score/hpo_max_score

        with open(output.part, "w") as w:
            w.write("set\tmethod\tparam\tmed_score\thpo_score\n")
            w.write(f"{wildcards.joinedname}\t{wildcards.method}\t{wildcards.param}\t{medicine_hitrate}\t{hpo_percent_score}\n")


rule aggregate_SC_cutoff_range:
    input:
        expand("work/edgelists/clustering/results_parts/{{joinedname}}_{cutoff}_SC.csv", cutoff = [round(x, 4) for x in np.arange(0, 1, 0.05, dtype=float)])
    output:
        "work/edgelists/clustering/results/{joinedname}_SC.csv"
    shell:
        """
        echo "set\tmethod\tparam\tmed_score\thpo_score" > {output}
        awk 'FNR > 1' {input} >> {output}
        """

rule aggregate_MCL_inflation_range:
    input:
        expand("work/edgelists/clustering/results_parts/{{joinedname}}_{inflation}_MCL.csv", inflation=[round(x, 4) for x in np.arange(1,5,0.2,dtype=float)])
    output:
        "work/edgelists/clustering/results/{joinedname}_MCL.csv"
    shell:
        """
        echo "set\tmethod\tparam\tmed_score\thpo_score" > {output}
        awk 'FNR > 1' {input} >> {output}
        """


rule aggregate_MCL_cutoff:
    input:
        expand("work/edgelists/clustering/results_parts/{{joinedname}}_{cutoff}_cutoffMCL.csv", cutoff=[round(x, 4) for x in np.arange(0,1,0.05,dtype=float)])
    output:
        "work/edgelists/clustering/results/{joinedname}_cutoffMCL.csv"
    shell:
        """
        echo "set\tmethod\tparam\tmed_score\thpo_score" > {output}
        awk 'FNR > 1' {input} >> {output}
        """