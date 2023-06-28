import numpy as np

include: "FormatEdgeList/FormatFilterNormalize.smk"
include: "Clustering/SpectralClustering.smk"
include: "Clustering/MCL.smk"
include: "Scoring/ClusterScore.smk"
include: "Plotting/Plot.smk"
include: "CandidateGenes/RankCandidateGenes.smk"
include: "Enrichment/Infer.smk"
include: "Enrichment/Enrichment.smk"

clusters = "work/edgelists/clustering/candidate_gene/enrichment/{joinedname}_{{cutoff}}_{{method}}_{{cluster}}_KEGG.csv".format(
            joinedname=config["name"])


wildcard_constraints:
    joinedname= "[a-zA-Z0-9\-]+",
    cutoff= "[0-9.]+",
    method= "[A-Za-z]+"
### Input all
rule all:
    input:
        "work/edgelists/plots/{joinedname}_SC.png". format(joinedname=config["name"]),
        "work/edgelists/plots/{joinedname}_MCL.png".format(joinedname=config["name"]),
        "work/edgelists/plots/{joinedname}_cutoffMCL.png".format(joinedname=config["name"]),
        expand(clusters,
            cutoff=[round(x, 4) for x in np.arange(0.75, 1, 0.05,dtype=float)],
            method=["SC", "cutoffMCL"],
            cluster=range(3)
        )