import pandas as pd

hpoid_df = pd.read_csv("data/micromedex/all_unique_terms.csv", sep="\t")
hpoid_df = hpoid_df[hpoid_df.HPO == hpoid_df.HPO]
hpoid_df = hpoid_df[hpoid_df.No_genes != hpoid_df.No_genes]

hpoid_dict = dict()
for i, row in hpoid_df.iterrows():
    term = row["Term"]
    hpoid = row["HPO"]
    hpoid_dict[term] = hpoid.replace(":","-")

cluster_files = [
    "anti_inflamation",
    "ssri_tricyclic_antidepressant",
    "statins"
    ]

HPO_lines = []
med_lines = []
skipped = ["tenoxicam", "lornoxicam"]
for clustername in cluster_files:
    cluster_file = f"data/selected_clusters/{clustername}.txt"
    med_in_cluster = [l.strip() for l in open(cluster_file, "r")]
    adr_count = dict()

    for med in med_in_cluster:
        if med not in skipped:
            adr_location = f"data/micromedex/{med.capitalize()}_adrs.txt"
            adrs = [l.strip() for l in open(adr_location, "r")]
            med_lines.append({
                "cluster": clustername,
                "drug": med.capitalize()
            })
            for adr in adrs:
                try:
                    hpoid = hpoid_dict[adr]
                    if hpoid in adr_count:
                        adr_count[hpoid] += 1
                    else:
                        adr_count[hpoid] = 1
                except KeyError:
                    pass

    for hpoid, count in adr_count.items():
        HPO_lines.append({
            "HPO": hpoid,
            "count": count,
            "cluster": clustername
        })

df_hpo = pd.DataFrame(HPO_lines)
summed_hpo = df_hpo.groupby("HPO")["count"].sum()
summed_hpo = summed_hpo.reset_index()
summed_hpo = summed_hpo.rename({"count": "count_total", "HPO": "HPO"}, axis=1)
df_hpo = df_hpo.merge(summed_hpo, on="HPO")
df_hpo["percent_count"] = df_hpo["count"]/df_hpo["count_total"]

df_hpo.to_csv("data/scoring/HPO-count-cluster.csv", sep="\t", index=False)
pd.DataFrame(med_lines).to_csv("data/scoring/med-cluster.csv", sep="\t", index=False)
