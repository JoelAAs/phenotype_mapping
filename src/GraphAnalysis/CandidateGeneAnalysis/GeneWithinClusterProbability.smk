import glob
import os
import bz2

project_name = config["project_name"]
terms = []
unique_genes = set()
for input_term in glob.glob(f"input/{project_name}/*csv"):
    term = os.path.basename(input_term).replace(".csv", "")
    terms.append(term)
    with open(input_term, "r") as f:
        genes = [l.strip() for l in f]
        config[term] = genes[1:]
        unique_genes.update(genes[1:])

config["terms"] = terms

checkpoint gene_in_cluster_probability_aggregation:
    input:
        #genes_in_paths = expand("work/{{project}}/neighborhood/{gene}_p_gene.csv.bz2", gene = unique_genes), problematic if its a joined project. hardcode and assume it exists
        term_genesets  = expand("input/{{project}}/{term}.csv", term = terms),
        cluster_file   = "work/{project}/clustering/SCnorm_{n_clusters}.csv"
    output:
        directory("work/{project}/cadidate_genes/probabilities_{n_clusters}/")
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
            probability_dict = {}
            terms = 0
            for node in cluster_dict[cluster]:
                project_input = ("HPO-pruned" if node[:3] == "HP-" or node[:5] == "ORPHA" else "full-drugbank")
                with open(f"input/{project_input}/{node}.csv", "r") as f:
                    genes = [l.strip() for l in f][1:]  # TODO: check if input has header
                    for gene in genes:
                        with bz2.open(f"work/{project_input}/neighborhood/{gene}_p_gene.csv.bz2", "r") as f:
                            header= True
                            for line in f:
                                if header:
                                    header = False
                                else:
                                    gene, p_gene_path = line.decode("utf-8").strip().split()
                                    if gene in probability_dict:
                                        probability_dict[gene] += p_gene_path
                                    else:
                                        probability_dict[gene] = p_gene_path
                terms += 1

            with open(f"work/{wildcards.project}/cadidate_genes/probabilities_{wildcards.n_clusters}/cluster_{cluster}.csv", "w") as w:
                w.write("gene\ty_probability\n")
                for gene in probability_dict:
                    w.write(f"{gene}\t{probability_dict[gene]/terms}\n")