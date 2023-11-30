include: "NetworkClusteringMethods/Clusters.smk"
include: "NetworkClusteringMethods/ClusterMetrics.smk"
include: "CandidateGeneAnalysis/GeneWithinClusterProbability.smk"
include: "NetworkPlotting/PlotGraphs.smk"

projects = [
#    "HPO-pruned",
    "full-drugbank",
    "MedAdr"
]
n_clusters = [
    2,
    3,
    4,
    5
]


## Rule
rule all:
    input:
        expand("work/{project}/plots/clustering/SC_{n_clusters}.png",
            project=projects, n_clusters=n_clusters),
        expand("work/{project}/plots/clustering/single_clusters_{n_clusters}/done.txt",
            project=projects, n_clusters=n_clusters),
        expand("work/{project}/plots/clustering/silhouette/silouette_{n_clusters}.png",
            project=projects, n_clusters=n_clusters),
        expand("work/{project}/candidate_genes/probabilities_{n}/",
            project=projects, n=[3,4]),
        expand("work/{project}/plots/candidate_genes/enrichment/enrichment_{n}/{m}/done.csv",
            project=projects, n=[3,4], m=["top", "diamond"]),
        expand("work/{project}/plots/candidate_genes/annotated/annotated_3/done.txt",
            project=projects)
