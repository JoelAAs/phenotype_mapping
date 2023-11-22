include: "NetworkClusteringMethods/Clusters.smk"
include: "NetworkPlotting/PlotGraphs.smk"

projects = [
    #"HPO-pruned",
    #"full-drugbank",
    "MedAdr"
]


## Rule
rule all:
    input:
        expand("work/{project}/clustering/plots/SC.png",
            project = projects),
        expand("work/{project}/clustering/plots/single_clusters/done.txt",
            project = projects)
