# Config
config["projects"] = [
    "HPO-pruned",
    "full-drugbank",
    "MedAdr"
]
n_clusters = [
    2,
    3,
    4
]

config["n_permut"] = 100
groups = ["Statin", "Antidepressant", "NSAID"]
config["n_permutation_combinations"] = 10000



include: "NetworkClusteringMethods/Clusters.smk"
include: "NetworkClusteringMethods/ClusterMetrics.smk"
include: "CandidateGeneAnalysis/GeneWithinClusterProbability.smk"
include: "NetworkPlotting/PlotGraphs.smk"
include: "CandidateGeneAnalysis/ClusterPermutation.smk"
include: "SeparationMeasure/Benchmarking.smk"
include: "CandidateGeneAnalysis/StandardScoreCandidateGenes.smk"
include: "CandidateGeneAnalysis/GeneProbabilityCDF.smk"

## Rule
rule all:
    input:
        expand("work/full-drugbank-benchmark/benchmarking/DIAMOnD_{group}.csv", group = groups),
        expand("work/full-drugbank-benchmark/group-quant/{group}_connected_components_enrichment/done.csv", group = groups)

