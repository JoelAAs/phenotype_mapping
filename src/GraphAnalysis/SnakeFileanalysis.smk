#### If multiproject project or singular
if "projects" not in config:
    config["projects"] = [config["project_name"]]

wildcard_constraints:
    n_clusters = "[0-9]+",
    cluster = "[0-9]+"

include: "NetworkClusteringMethods/Clusters.smk"
include: "NetworkClusteringMethods/ClusterMetrics.smk"
include: "CandidateGeneAnalysis/GeneWithinClusterProbability.smk"
include: "NetworkPlotting/PlotGraphs.smk"
include: "CandidateGeneAnalysis/ClusterPermutation.smk"
include: "CandidateGeneAnalysis/StandardScoreCandidateGenes.smk"
include: "CandidateGeneAnalysis/GeneProbabilityCDF.smk"



## Rule
rule all:
    input:
        expand("work/{project}/group-quant_{n_clusters}/{cluster}_connected_components_enrichment/done.csv",
            project=config["project_name"],
            n_clusters=config["n_clusters"],
            cluster=range(config["n_clusters"])
        )

