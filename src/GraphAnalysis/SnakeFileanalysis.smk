include: "NetworkClusteringMethods/Clusters.smk"
include: "NetworkClusteringMethods/ClusterMetrics.smk"
include: "CandidateGeneAnalysis/GeneWithinClusterProbability.smk"
include: "NetworkPlotting/PlotGraphs.smk"

projects = [
    "HPO-pruned",
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
        expand("work/{project}/clustering/plots/SC_{n_clusters}.png",
            project=projects, n_clusters=n_clusters),
        expand("work/{project}/clustering/plots/single_clusters_{n_clusters}/done.txt",
            project=projects, n_clusters=n_clusters),
        expand("work/{project}/clustering/plots/silouette_{n_clusters}.png",
            project=projects, n_clusters=n_clusters),
        expand("work/{project}/candidate_genes/probabilities_{n}/",
            project=projects, n=[3,4]),
        expand("work/{project}/candidate_genes/enrichment_{n}/{m}/done.csv",
            project=projects, n=[3,4], m=["top", "diamond"]),
        expand("work/{project}/candidate_genes/annotated_3/done.txt",
            project=projects)
