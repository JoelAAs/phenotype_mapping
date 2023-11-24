include: "NetworkClusteringMethods/Clusters.smk"
include: "NetworkPlotting/PlotGraphs.smk"
include: "NetworkClusteringMethods/ClusterMetrics.smk"
include: "CandidateGeneAnalysis/GeneWithinClusterProbability.smk"
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
        "work/MedAdr/cadidate_genes/probabilities_4/",
        "work/MedAdr/cadidate_genes/probabilities_3/",
        "work/HPO-pruned/cadidate_genes/probabilities_4/",
        "work/full-drugbank/cadidate_genes/probabilities_3/"