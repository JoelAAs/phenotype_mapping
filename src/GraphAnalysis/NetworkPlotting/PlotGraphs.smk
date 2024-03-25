from plot_graph import graph_plot, cluster_graph_plot, plot_silhouette
import pandas as pd
import os
import glob

def get_clusteredge_output(wc):
    cluster_edge_ck = checkpoints.ClusterEdges.get(**wc).output[0]
    clusters, = glob_wildcards(os.path.join(cluster_edge_ck, "SC_edges_{cluster}.csv"))
    return expand(os.path.join(cluster_edge_ck, "SC_edges_{cluster}.csv"), cluster=clusters)


rule SCPlot:
    input:
        graph = "work/{project}/term_to_term_probability_matrix.csv",
        clusters = "work/{project}/clustering/SCnorm_{n_cluster}.csv"
    output:
        figure = "work/{project}/plots/clustering/SC_{n_cluster}.png"
    run:
        graph_plot(
            input.graph,
            input.clusters,
            wildcards.project,
            "Spectral Clustering",
            output.figure)

checkpoint ClusterEdges:
    input:
        graph="work/{project}/term_to_term_probability_matrix.csv",
        clusters="work/{project}/clustering/SCnorm_{n_clusters}.csv"
    output:
        directory("work/{project}/clustering/clustered_edges_{n_clusters}")
    run:
        os.mkdir(output[0])

        cluster_dict = {}
        with open(input.clusters) as f:
            lines = f.readlines()
            for line in lines[1:]:
                node, cluster = line.strip().split("\t")
                if cluster in cluster_dict:
                    cluster_dict[cluster].append(node)
                else:
                    cluster_dict[cluster] = [node]

        edge_df = pd.read_csv(input.graph, sep = "\t")
        for cluster in cluster_dict:
            current_nodes = cluster_dict[cluster]

            cluster_df = edge_df[edge_df["query"].isin(current_nodes)]
            cluster_df = cluster_df[cluster_df["neighborhood"].isin(current_nodes)]
            cluster_df.to_csv(
                f"work/{wildcards.project}/clustering/clustered_edges_{wildcards.n_clusters}/SC_edges_{cluster}.csv",
                sep="\t",
                index=None,
            )

rule ClusterPlot:
    input:
        get_clusteredge_output
    output:
        check = "work/{project}/plots/clustering/single_clusters_{n_clusters}/done.txt"
    run:
        with open(output.check, "w") as w:
            w.write("check")
        for i, cluster_edge_file in enumerate(input):
            cluster_graph_plot(
                cluster_edge_file,
                wildcards.project,
                "Spectral clustering",
                f"work/{wildcards.project}/plots/clustering/single_clusters_{wildcards.n_clusters}/cluster_{i}.png"
            )

rule SilhouettePlot:
    input:
        silhouette_file = "work/{project}/clustering/metrics/silouette_{n_clusters}.csv"
    output:
        figure_location = "work/{project}/plots/clustering/silhouette/silouette_{n_clusters}.png"
    run:
        plot_silhouette(
            input.silhouette_file,
            output.figure_location
        )

rule ScatterPlot:
    input:
        probability_annotatied_csv = "work/{project}/candidate_genes/annotated_{n_clusters}/annotated_{cluster}.csv"
    output:
        figure = "work/{project}/plots/candidate_genes/annotated/annotated_{n_clusters}/annotated_{cluster}.png"
    shell:
        """
        Rscript src/GraphAnalysis/NetworkPlotting/PlotClusterGeneProbability.R {input} {output}
        """



rule PlotEnrichment:
    input:
        enrichments = "work/{project}/candidate_genes/enrichment_{n_clusters}/{method}/enrichment_{cluster}.csv"
    output:
        figure = "work/{project}/plots/candidate_genes/enrichment/enrichment_{n_clusters}/{method}/enrichment_{cluster}.png"
    shell:
        """
        Rscript src/GraphAnalysis/NetworkPlotting/Enrichmentplot.R {input} {output}
        """

rule PlotEnrichment_LCC:
    input:
        enrichments = "work/{project}/group-quant_{n_clusters}/{cluster}_connected_components_enrichment/Component_{component}_KEGG.csv"
    output:
        figure = "work/{project}/group-quant_{n_clusters}/{cluster}_connected_components_enrichment/Component_{component}_KEGG.png"
    shell:
        """
        Rscript src/GraphAnalysis/NetworkPlotting/Enrichmentplot.R {input} {output}
        """



rule PlotPermutations:
    input:
        drug_permut = "work/full-drugbank/clustering/metrics/sensitivity_inbew.csv",
        drug_permut_norm= "work/full-drugbank/clustering/metrics/sensitivity_norm.csv",
        sab_drug_permut = "work/full-drugbank-benchmark/clustering/metrics/sensitivity_sab.csv",
        zscore_drug_permut = "work/full-drugbank-benchmark/clustering/metrics/sensitivity_zscore.csv",
        adr_drug_permut = "work/MedAdr/clustering/metrics/sensitivity_inbew.csv",
        adr_drug_permut_norm= "work/MedAdr/clustering/metrics/sensitivity_norm.csv",
        sab_adr_drug_permut = "work/MedAdr-benchmark/clustering/metrics/sensitivity_sab.csv",
        zscore_adr_drug_permut = "work/MedAdr-benchmark/clustering/metrics/sensitivity_zscore.csv",
        adr_permut= "work/MedAdr/clustering/metrics/hpo_sensitivity_inbew.csv",
        adr_permut_norm= "work/MedAdr/clustering/metrics/hpo_sensitivity_norm.csv",
        sab_adr_permut = "work/MedAdr-benchmark/clustering/metrics/hpo_sensitivity_sab.csv",
        zscore_adr_permut= "work/MedAdr-benchmark/clustering/metrics/hpo_sensitivity_zscore.csv"
    output:
        "work/meta_plots/Sensitivity.png"
    shell:
        """
        Rscript src/GraphAnalysis/NetworkPlotting/NodeSensitivityPermutationPlot.R {input} {output}
        """

