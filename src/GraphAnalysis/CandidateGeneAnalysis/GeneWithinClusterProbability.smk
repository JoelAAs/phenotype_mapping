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


def get_expected_enrichments(wc):
    cluster_prob_ck = checkpoints.gene_in_cluster_probability_aggregation.get(**wc).output[0]
    clusters, = glob_wildcards(os.path.join(cluster_prob_ck, "cluster_{cluster}.csv"))
    output_path = f"work/{wc.project}/plots/candidate_genes/enrichment/enrichment_{wc.n_clusters}/{wc.method}/enrichment_{{cluster}}.png"
    return expand(output_path, cluster = clusters)




rule gene_in_cluster_probability_aggregation:
    input:
        genes_in_paths = get_gene_input(config["projects"]),
        term_genesets  = get_term_genesets(config["projects"]),
        cluster_file   = "work/{project}/clustering/SCnorm_{n_clusters}.csv"
    output:
        prob_dir = expand("work/{{project}}/candidate_genes/probabilities_{{n_clusters}}/cluster_{cluster}.csv",
            cluster=range(config["n_clusters"]))
    run:
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
            probability_dict = dict()
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


rule per_cluster_unique_input:
    input:
        term_genesets  = get_term_genesets(config["projects"]),
        cluster_file   = "work/{project}/clustering/SCnorm_{n_clusters}.csv"
    output:
        unique_dir = expand("work/{{project}}/candidate_genes/unique_input_{{n_clusters}}/unique_input_{cluster}.csv",
            cluster=range(config["n_clusters"]))
    run:
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
            unique_input_genes = set()
            for node in cluster_dict[cluster]:
                for project in config["projects"]:
                    if os.path.exists(f"input/{project}/{node}.csv"):
                        project_input = project
                        break

                with open(f"input/{project_input}/{node}.csv", "r") as f:
                    genes = [l.strip() for l in f][1:]
                    for gene in genes:
                        unique_input_genes.add(gene)

            with open(f"work/{wildcards.project}/candidate_genes/unique_input_{wildcards.n_clusters}" + f"/unique_input_{cluster}.csv", "w") as w:
                w.write("gene\n")
                for gene in unique_input_genes:
                    w.write(f"{gene}\n")

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


rule get_norm_term_to_term_matrix:
    input:
        hpo = "work/HPO-pruned/term_to_term_probability_matrix.csv",
        drug = "work/full-drugbank/term_to_term_probability_matrix.csv"
    output:
        hpo_norm = "work/HPO-pruned/term_to_term_probability_matrix_norm.csv",
        drug_norm = "work/full-drugbank/term_to_term_probability_matrix_norm.csv"
    run:
        hpo_df = pd.read_csv(input.hpo, sep = "\t")
        hpo_df["probability"] = hpo_df["probability"]/max(hpo_df["probability"])
        hpo_df.to_csv(output.hpo_norm, sep="\t", index=None)

        drug_df = pd.read_csv(input.drug,sep="\t")
        drug_df["probability"] = drug_df["probability"] / max(drug_df["probability"])
        drug_df.to_csv(output.drug_norm,sep="\t",index=None)