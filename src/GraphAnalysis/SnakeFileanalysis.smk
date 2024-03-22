include: "NetworkClusteringMethods/Clusters.smk"
include: "NetworkClusteringMethods/ClusterMetrics.smk"
include: "CandidateGeneAnalysis/GeneWithinClusterProbability.smk"
include: "NetworkPlotting/PlotGraphs.smk"
include: "CandidateGeneAnalysis/ClusterPermutation.smk"
include: "CandidateGeneAnalysis/StandardScoreCandidateGenes.smk"
include: "CandidateGeneAnalysis/GeneProbabilityCDF.smk"

#### If multiproject project or singular
if "projects" not in config:
    config["projects"] = [config["project_name"]]


## Rule
rule all:
    input:
        expand("work/{project}/candidate_genes/probabilities_{n_clusters}", project=config["project_name"], n_clusters=config["clusters"]),
        expand("work/{project}/group-quant_{n_clusters}/{cluster}_quant.csv",
            project=config["project_name"],
            n_clusters=config["clusters"],
            cluster=range(config["clusters"])
        )


            #expand("work/full-drugbank-benchmark/group-quant/{group}_connected_components_enrichment/done.csv", group = groups)

