import glob
import os
import bz2
import pandas as pd
import numpy as np
# from src.GraphAnalysis.SeparationMeasure.DIAMOnD import *
from ClusterGeneProbabilities import get_centrality_of_cluster_genes, get_occurrence_of_genes_in_terms

# Config

for project in config["projects"]:
    terms = glob.glob(f"input/{project}/*.csv")
    terms = [t.split("/")[2].replace(".csv", "") for t in terms]
    config[project] = terms


# input functions
def get_gene_input(projects):
    expected_gene_input = []
    for project in projects:
        for term in config[project]:
            with open(f"input/{project}/{term}.csv", "r") as f:
                genes = [g.strip() for g in f][1:]
                for gene in genes:
                    expected_gene_input.append(f"work/{project}/neighborhood/{gene}_p_gene.csv.bz2")
    return list(set(expected_gene_input))

def get_term_genesets(projects):
    expected_term_genests = []
    for project in projects:
        for term in config[project]:
            expected_term_genests.append(f"input/{project}/{term}.csv")
    return list(set(expected_term_genests))

def get_cluster_gene_probabilities(wc):
    cluster_prob_ck = checkpoints.gene_in_cluster_probability_aggregation.get(**wc).output[0]
    clusters, = glob_wildcards(os.path.join(cluster_prob_ck, "cluster_{cluster}.csv"))
    return expand(os.path.join(cluster_edge_ck, "cluster_{cluster}.csv"), cluster=clusters)

def get_expected_enrichments(wc):
    cluster_prob_ck = checkpoints.gene_in_cluster_probability_aggregation.get(**wc).output[0]
    clusters, = glob_wildcards(os.path.join(cluster_prob_ck, "cluster_{cluster}.csv"))
    output_path = f"work/{wc.project}/plots/candidate_genes/enrichment/enrichment_{wc.n_clusters}/{wc.method}/enrichment_{{cluster}}.png"
    return expand(output_path, cluster = clusters)

def get_expected_centrality(wc):
    cluster_prob_ck = checkpoints.gene_in_cluster_probability_aggregation.get(**wc).output[0]
    clusters, = glob_wildcards(os.path.join(cluster_prob_ck,"cluster_{cluster}.csv"))
    output_path = f"work/{wc.project}/plots/candidate_genes/annotated/annotated_{wc.n_clusters}/annotated_{{cluster}}.png"
    return expand(output_path,cluster=clusters)


checkpoint gene_in_cluster_probability_aggregation:
    input:
        genes_in_paths = get_gene_input(config["projects"]),
        term_genesets  = get_term_genesets(config["projects"]),
        cluster_file   = "work/{project}/clustering/SCnorm_{n_clusters}.csv"
    output:
        directory("work/{project}/candidate_genes/probabilities_{n_clusters}/")
    run:
        os.mkdir(output[0])
        cluster_dict = {}
        with open(input.cluster_file,"r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                node, cluster = line.strip().split("\t")
                if cluster in cluster_dict:
                    cluster_dict[cluster].append(node)
                else:
                    cluster_dict[cluster] = [node]

        for cluster in cluster_dict:
            probability_dict = {}
            for node in cluster_dict[cluster]:
                print(f"Node: {node} for cluster {cluster}")
                for project in config["projects"]:
                    if os.path.exists(f"input/{project}/{node}.csv"):
                        project_input = project
                        break

                with open(f"input/{project_input}/{node}.csv", "r") as f:
                    genes = [l.strip() for l in f][1:]
                    for gene in genes:
                        with bz2.open(f"work/{project_input}/neighborhood/{gene}_p_gene.csv.bz2", "r") as f:
                            header= True
                            for line in f:
                                if header:
                                    header = False
                                else:
                                    gene, p_gene_path = line.decode("utf-8").strip().split()
                                    p_gene_path = float(p_gene_path)/len(genes)
                                    if gene in probability_dict:
                                        probability_dict[gene] += p_gene_path
                                    else:
                                        probability_dict[gene] = p_gene_path

            with open(f"work/{wildcards.project}/candidate_genes/probabilities_{wildcards.n_clusters}/cluster_{cluster}.csv", "w") as w:
                w.write("gene\ty_probability\n")
                for gene in probability_dict:
                    w.write(f"{gene}\t{probability_dict[gene]/len(cluster_dict[cluster])}\n")


rule get_stringid:
    params:
        centrality_max = 0.3
    input:
        cluster_genes = "work/{project}/candidate_genes/annotated_{n_clusters}/annotated_{cluster}.csv"
    output:
        top = "work/{project}/candidate_genes/enrichment_{n_clusters}/diamond/input/string_top_{cluster}.csv"
    run:
        cluster_prob = pd.read_csv(input.cluster_genes, sep="\t")

        cluster_prob["centrality_norm"] = cluster_prob["centrality"]/cluster_prob["centrality"].max()


        cluster_prob["logprob"] = np.log10(cluster_prob["y_probability"])  # approx normal dist in logspace
        mu = cluster_prob["logprob"].mean()
        std = cluster_prob["logprob"].std()

        threshold = mu + std * 2
        top_df = cluster_prob[cluster_prob["logprob"] > threshold]
        cluster_prob = cluster_prob[cluster_prob["centrality_norm"] < params.centrality_max]
        top_df.to_csv(output.top, index=False, sep = "\t")



# rule Diamond:
#     params:
#         n_genes= 50
#     input:
#         cluster_genes = "work/{project}/candidate_genes/enrichment_{n_clusters}/diamond/input/string_top_{cluster}.csv",
#         ppi_cluster = "data/9606.protein.links_above_700.v11.5.txt",
#     output:
#         inferred_genes = "work/{project}/candidate_genes/enrichment_{n_clusters}/diamond/enriched/top_stringid_{cluster}.csv"
#     run:
#         input_list = ["", input.ppi_cluster, input.cluster_genes, params.n_genes, 1, output.inferred_genes]
#         network_edgelist_file, seeds_file, max_number_of_added_nodes, alpha, outfile_name = check_input_style(input_list)
#
#         G_original, seed_genes = read_input(network_edgelist_file, seeds_file)
#
#         added_nodes = DIAMOnD(
#             G_original,
#             seed_genes,
#             max_number_of_added_nodes, alpha,
#             outfile=outfile_name)
#
#
# rule translate_diamond:
#     input:
#         diamond_top = "work/{project}/candidate_genes/enrichment_{n_clusters}/diamond/enriched/top_stringid_{cluster}.csv",
#         seed_tops = "work/{project}/candidate_genes/enrichment_{n_clusters}/diamond/input/string_top_{cluster}.csv",
#         entrez= "data/ncbi/entrez.csv",
#         string_id="data/stringdb/9606.protein.info.v11.5.txt",
#     output:
#         string_top = "work/{project}/candidate_genes/enrichment_{n_clusters}/diamond/top_{cluster}.csv",
#         translated_preffered= "work/{project}/candidate_genes/enrichment_{n_clusters}/top/top_gene_name{cluster}.csv"
#     run:
#         diamond_top = pd.read_csv(input.diamond_top, sep = "\t")
#         del diamond_top["#rank"]
#         del diamond_top["p_hyper"]
#         diamond_top = diamond_top.rename({"DIAMOnD_node": "gene"}, axis = 1)
#         seed_top = pd.read_csv(input.seed_tops, sep = "\t")
#         del seed_top["y_probability"]
#         del seed_top["logprob"]
#         all_top = pd.concat([diamond_top, seed_top])
#
#         entrez_df = pd.read_csv(input.entrez,sep="\t")
#         string_df = pd.read_csv(input.string_id,sep="\t")
#
#         string_top = diamond_top.merge(
#             string_df,
#             left_on="gene",
#             right_on="string_protein_id",
#             how="left")
#
#         entrez_top = string_top.merge(
#             entrez_df,
#             left_on="preferred_name",
#             right_on="gene_name",
#             how="left")
#
#         entrez_top["entrez"].to_csv(output.string_top,sep="\t",index=False)

rule probability_cutoff_and_entrez:
    params:
        centrality_max = 0.3,
        exclude_starting_terms = True,
        top_n = 200
    input:
        entrez = "data/ncbi/entrez.csv",
        string_id= "data/stringdb/9606.protein.info.v11.5.txt",
        cluster_probability= "work/{project}/candidate_genes/annotated_{n_clusters}/annotated_{cluster}.csv"
    output:
        translated = "work/{project}/candidate_genes/enrichment_{n_clusters}/top/top_{cluster}.csv",
        translated_preffered = "work/{project}/candidate_genes/enrichment_{n_clusters}/top/top_gene_name{cluster}.csv",
    run:
        cluster_prob = pd.read_csv(input.cluster_probability, sep = "\t")
        #cluster_prob["centrality_norm"] = cluster_prob["centrality"]/cluster_prob["centrality"].max()


        # cluster_prob["logprob"] = np.log10(cluster_prob["y_probability"]) # approx normal dist in logspace

        #cluster_prob = cluster_prob[cluster_prob["logprob"] > threshold]
        #cluster_prob = cluster_prob[cluster_prob["centrality_norm"] < params.centrality_max]

        if params.exclude_starting_terms:
            cluster_prob = cluster_prob[cluster_prob["count"] == 0]

        cluster_prob = cluster_prob.sort_values(by="y_probability", ascending=False)[:params.top_n]

        entrez_df = pd.read_csv(input.entrez, sep="\t")
        string_df = pd.read_csv(input.string_id, sep = "\t")

        string_top = cluster_prob.merge(
            string_df,
            left_on = "gene",
            right_on="string_protein_id",
            how="left")

        entrez_top = string_top.merge(
            entrez_df,
            left_on="preferred_name",
            right_on="gene_name",
            how="left")

        entrez_top[["gene_name", "string_protein_id", "entrez"]].to_csv(output.translated_preffered, sep = "\t", index=False)
        entrez_top["entrez"].to_csv(output.translated, sep="\t", index=False)

rule enrichment_analysis:
    input:
        cluster_probability = "work/{project}/candidate_genes/enrichment_{n_clusters}/{method}/top_{cluster}.csv"
    output:
        enrichment_results = "work/{project}/candidate_genes/enrichment_{n_clusters}/{method}/enrichment_{cluster}.csv"
    shell:
        """
        Rscript src/GraphAnalysis/CandidateGeneAnalysis/Enrichment.R {input} {output}
        """


rule enrichment_done:
    input:
        get_expected_enrichments
    output:
        "work/{project}/plots/candidate_genes/enrichment/enrichment_{n_clusters}/{method}/done.csv"
    shell:
        """
        touch {output}
        """

rule get_centrality_and_frequence:
    input:
        probabilities_df = "work/{project}/candidate_genes/probabilities_{n_clusters}/cluster_{cluster}.csv",
        cluster_file = "work/{project}/clustering/SCnorm_{n_clusters}.csv",
        ppi_file = "data/9606.protein.links.v11.5.txt",
    output:
        probability_annotatied_csv = "work/{project}/candidate_genes/annotated_{n_clusters}/annotated_{cluster}.csv"
    run:
        term_frequence_df, n_terms = get_occurrence_of_genes_in_terms(
            input.cluster_file,
            wildcards.cluster,
            config["projects"]
        )
        centrality_df = get_centrality_of_cluster_genes(
            input.ppi_file,
            input.probabilities_df
        )
        centrality_df["n_terms"] = n_terms

        annotated_df = centrality_df.merge(
            term_frequence_df,
            on="gene",
            how="left"
        )
        annotated_df = annotated_df.fillna(value=0)
        annotated_df.to_csv(
            output.probability_annotatied_csv, sep="\t", index=False
        )

rule get_centrality_and_frequnece_done:
    input:
        get_expected_centrality
    output:
        "work/{project}/plots/candidate_genes/annotated/annotated_{n_clusters}/done.txt"
    shell:
        """
        touch {output}
        """