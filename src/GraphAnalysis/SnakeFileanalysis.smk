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
            project=projects, n_clusters=n_clusters)
        #"work/MedAdr/candidate_genes/probabilities_4/",
        #"work/MedAdr/candidate_genes/probabilities_3/",
        #"work/MedAdr/candidate_genes/enrichment_4/top/done.csv",
        #"work/MedAdr/candidate_genes/enrichment_4/diamond/done.csv"
