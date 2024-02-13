import bz2
import numpy as np

group_dict = dict()
with open("data/Node_groups.csv", "r") as f:
    for l in f.readlines()[1:]:
        drug, group = l.strip().split("\t")
        if group in group_dict:
            group_dict[group].append(drug)
        else:
            group_dict[group] = [drug]

for group, drugs in group_dict.items():
    config[group] = drugs


def get_genes_in_permutset(wc):
    with open(f"input/{wc.project}/{wc.term}.csv", "r") as f:
        genes = [l.strip() for l in f.readlines()[1:]]

    return  genes


rule probabilities_of_gene_in_term:
    input:
        term = "input/{project}/{term}.csv",
        gene_probabilites = lambda wc: expand("work/{{project}}/neighborhood/{gene}_p_gene.csv.bz2", gene = get_genes_in_permutset(wc))
    output:
        term_prob = "work/{project}/candidate_genes/term_probabilities/{term}_probabilities.csv"
    run:
        n_genes = 0
        prob_dict = dict()
        for gene_prob in input.gene_probabilites:
            with bz2.open(gene_prob, "r") as f:
                n_genes += 1
                for l in f.readlines()[1:]:
                    gene, prob = l.decode("utf-8").strip().split("\t")
                    prob = float(prob)
                    if gene in prob_dict:
                        prob_dict[gene] += prob
                    else:
                        prob_dict[gene] = prob

        with open(output.term_prob,"w") as w:
            w.write(f"gene\tprobability\n")
            for gene, prob in prob_dict.items():
                w.write(f"{gene}\t{prob/n_genes}\n")


rule get_mean_and_variation_per_term:
    input:
        term_prob = expand(
            "work/{{project}}/candidate_genes/term_probabilities/{{term}}_set_{n}_probabilities.csv",
            n = range(config["n_permut"]))
    output:
        permut_mean_var = "work/{project}/candidate_genes/term_probabilities_metrics/{term}_mean_var.csv"
    run:
        permut_prob_dict = dict()
        for i, permut_probabilities in enumerate(input.term_prob):
            with open(permut_probabilities) as f:
                for l in f.readlines()[1:]:
                    gene, prob = l.strip().split("\t")
                    if gene not in permut_prob_dict:
                        permut_prob_dict[gene] = np.zeros(config["n_permut"])
                    permut_prob_dict[gene][i] = float(prob)

        with open(output.permut_mean_var, "w") as w:
            w.write("gene\tmean\tstd\n")
            for gene, probs in permut_prob_dict.items():
                mean = probs.mean()
                var = probs.var()
                w.write(f"{gene}\t{mean}\t{var}\n")

rule get_mean_and_var_cluster:
    input:
        term_mean_var = lambda wc: expand(
            "work/{{project}}/candidate_genes/term_probabilities_metrics/{term}_mean_var.csv",
            term = config[wc.group])
    output:
        cluster_prob = "work/{project}/candidate_genes/term_probabilities_groups/{group}_mean_var.csv"
    run:
        permut_prob_dict = dict()
        for i, permut_probabilities in enumerate(input.term_mean_var):
            with open(permut_probabilities) as f:
                for l in f.readlines()[1:]:
                    gene, mean, variance = l.strip().split("\t")
                    if gene not in permut_prob_dict:
                        permut_prob_dict[gene] = {
                            "mean": np.zeros(config["n_permut"]),
                            "variance": np.zeros(config["n_permut"])
                        }
                    permut_prob_dict[gene]["mean"][i] = float(mean)
                    permut_prob_dict[gene]["variance"][i] = float(variance)

        with open(output.cluster_prob, "w") as w:
            w.write("gene\tmean\tstd\n")
            for gene, probs_variance in permut_prob_dict.items():
                mean = probs_variance["mean"].mean() # Assumed normal distribution
                var = sum(probs_variance["variance"])*len(probs_variance["variance"])**-2
                w.write(f"{gene}\t{mean}\t{var}\n")



rule calculate_z_scores:
    input:
