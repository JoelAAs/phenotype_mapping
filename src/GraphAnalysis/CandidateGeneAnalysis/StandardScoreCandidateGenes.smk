import bz2
import readline

import numpy as np
import random


def get_genes_in_permutset(wc):
    with open(f"input/{wc.project}/{wc.term}.csv", "r") as f:
        genes = [l.strip() for l in f.readlines()[1:]]

    return  genes

def get_terms_in_cluster(wc):
    group_dict = dict()
    with open(f"work/{config['project_name']}/clustering/SCnorm_{config['n_clusters']}.csv", "r") as f:
        for l in f.readlines()[1:]:
            node, cluster = l.strip().split("\t")
            if cluster not in group_dict:
                group_dict[cluster] = []
            group_dict[cluster].append(node)

    return expand(
        "work/{permutation_project}/candidate_genes/term_probabilities/{term}_set_{n}_probabilities.csv",
        permutation_project=config["permutation_folder"],
        term=group_dict[wc.cluster],
        n=config["permutation_N"]
    )

rule get_permuted_cluster_probabilities:
    params:
        number_of_combinations = config["n_permutation_combinations"]
    input:
        clusters = "work/{project}/clustering/SCnorm_{n_clusters}.csv".format(
            project=config["project_name"],
            n_clusters=config["n_clusters"]
        ),
        term_probabilities = lambda wc: get_terms_in_cluster(wc)
    output:
        cluster_prob = "work/{permutation_project}/candidate_genes/cluster_gene_permutations_{n_clusters}/{cluster}_probability.csv"
    run:
        group_dict = dict()
        with open(input.clusters, "r") as f:
            for l in f.readlines()[1:]:
                node, cluster = l.strip().split("\t")
                if cluster not in group_dict:
                    group_dict[cluster] = []
                group_dict[cluster].append(node)

        group_terms = group_dict[wildcards.cluster]

        permut_prob_dict = dict()
        permutations_sets = set()
        while len(permutations_sets) < params.number_of_combinations:
            permutations_sets.add(tuple(random.sample(range(config["n_permut"]), k=len(group_terms))))
        k = 0
        for permut_set in permutations_sets:
            print(f"{k} combination of {params.number_of_combinations}")
            k += 1
            for j, (term, i) in enumerate(zip(group_terms, list(permut_set))):
                term_permut_file = f"work/full-drugbank-benchmark/candidate_genes/term_probabilities/{term}_set_{i}_probabilities.csv"
                with open(term_permut_file, "r") as f:
                    for l in f.readlines()[1:]:
                        gene, prob = l.strip().split("\t")
                        prob = float(prob)

                        if gene not in permut_prob_dict:
                            permut_prob_dict[gene] = np.zeros(params.number_of_combinations)

                        permut_prob_dict[gene][j] += prob/len(group_terms)

        with open(output.cluster_prob, "w") as w:
            genes = list(permut_prob_dict.keys())
            w.write("gene\tprob\tn\n")
            for gene in genes:
                probs = permut_prob_dict[gene]
                for i, prob in enumerate(probs):
                    w.write(f"{gene}\t{prob}\t{i}\n")






