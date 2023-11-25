import glob
import os
import bz2
import numpy as np
import pandas as pd

# Config
# if type(config["projects"]) != list():
#     config["projects"] = [config["projects"], ]


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
    output_path = f"work/{wc.project}/clustering/plots/enrichment_{wc.n_clusters}/{wc.method}/enrichment_{{cluster}}.png"
    return expand(output_path, cluster = clusters)

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
            terms = 0
            for node in cluster_dict[cluster]:
                print(f"Node: {node} for cluster {cluster}")
                project_input = ("HPO-pruned" if node[:3] == "HP-" or node[:5] == "ORPHA" else "full-drugbank")
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
                                    p_gene_path = float(p_gene_path)
                                    if gene in probability_dict:
                                        probability_dict[gene] += p_gene_path
                                    else:
                                        probability_dict[gene] = p_gene_path
                terms += 1

            with open(f"work/{wildcards.project}/candidate_genes/probabilities_{wildcards.n_clusters}/cluster_{cluster}.csv", "w") as w:
                w.write("gene\ty_probability\n")
                for gene in probability_dict:
                    w.write(f"{gene}\t{probability_dict[gene]/terms}\n")

rule probability_cutof_and_entrez:
    input:
        entrez = "data/ncbi/entrez.csv",
        string_id= "data/stringdb/9606.protein.info.v11.5.txt",
        cluster_probability= "work/{project}/candidate_genes/probabilities_{n_clusters}/cluster_{cluster}.csv"
    output:
        translated = "work/{project}/candidate_genes/enrichment_{n_clusters}/top/entrez_top_{cluster}.csv"
    run:
        cluster_prob = pd.read_csv(input.cluster_probability, sep = "\t")
        cluster_prob["logprob"] = np.log10(cluster_prob["y_probability"]) # approx normal dist in logspace
        mu = cluster_prob["logprob"].mean()
        std = cluster_prob["logprob"].std()

        threshold = mu + std*2
        top_df = cluster_prob[cluster_prob["logprob"] > threshold]

        entrez_df = pd.read_csv(input.entrez, sep="\t")
        string_df = pd.read_csv(input.string_id, sep = "\t")

        string_top = top_df.merge(
            string_df,
            left_on = "gene",
            right_on="string_protein_id",
            how="left")

        entrez_top = string_top.merge(
            entrez_df,
            left_on="preferred_name",
            right_on="gene_name",
            how="left")

        entrez_top["entrez"].to_csv(output.translated, index=False)


rule Diamond:
    params:
        ppi_cutoff=0.7,
    input:
        cluster_genes = "work/{project}/candidate_genes/enrichment_{n_clusters}/top/entrez_top_{cluster}.csv",
        ppi_cluster = "data/9606.protein.links.v11.5.txt",
    output:
        inferred_genes = "work/{project}/candidate_genes/enrichment_{n_clusters}/diamond/entrez_top_{cluster}.csv"
    singularity:
        "envs/modifier.simg"
    shell:
        """
        Rscript src/GraphAnalysis/CandidateGeneAnalysis/Diamond.R {input.cluster_genes} {output.inferred_genes} \
            {input.ppi_cluster} {params.ppi_cutoff}
        """

rule enrichment_analysis:
    input:
        cluster_probability = "work/{project}/candidate_genes/enrichment_{n_clusters}/{method}/entrez_top_{cluster}.csv"
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
        "work/{project}/candidate_genes/enrichment_{n_clusters}/{method}/done.csv"
    shell:
        """
        touch {output}
        """