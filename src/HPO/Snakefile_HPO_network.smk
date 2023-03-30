import sys
import pandas as pd
from GetHPOGenes import GetHPOGenes
from testHPO_paralleled import GetGeneDistanceAtN
from itertools import combinations
import re
import os
import bz2

def  get_all_hpos(file_csv):
    df = pd.read_csv(file_csv, sep = "\t")
    hpos = set(list(df.hpoid))
    hpos = [h for h in hpos if h == h]
    current_hpo_combinations = list(combinations(hpos, 2))
    sorted_hpo_combinations = [sorted(comb) for comb in current_hpo_combinations]
    return sorted_hpo_combinations, hpos

hpo_combinations, all_hpos = get_all_hpos("data/MatchedHPO.csv")
combination_input = ["work/HPO/distances/{}_{}_all.csv.bz2".format(comb[0], comb[1]) for comb in hpo_combinations]

wildcard_constraints:
    hpo="[A-Za-z0-9.: \-()]+"





### Input functions

def get_all_genes_at_depth(wc):

    all_needed_files_exists = all(
        map(
            os.path.exists,
        [
            "work/HPO/gene_path",
            "work/HPO/pairs/name_conversion.csv",
            "work/HPO/pairs/distance_pairs.csv"]))
    if not all_needed_files_exists:
        re_eval = checkpoints.get_possible_paths_from_genes.get()
    else:
        df_distance_pairs = pd.read_csv("work/HPO/pairs/distance_pairs.csv", sep = "\t", header=None)
        expected_distance = []
        for i, gene_a, gene_b, depth in df_distance_pairs.itertuples():
            expected_distance.append(
                f"work/HPO/gene_path/{gene_a}_at_{depth}_shortest.csv.bz2"
            )
            expected_distance.append(
                f"work/HPO/gene_path/{gene_b}_at_{depth}_shortest.csv.bz2"
            )

        return expected_distance

# rule get_scores:
#     input:
#         distances = combination_input
#     output:
#         summed_scores = "work/chembl/distances/summed_scores.csv"
#     run:
#         summed_scores = []
#         distance_files = input.distances
#         nr_files = len(distance_files)
#         for i, distance_file in enumerate(input.distances):
#             print(f"{round(i/nr_files,3)*100} % done:", end="\r")
#             df_scores = pd.read_csv(distance_file, sep = "\t")
#             df_scores["score_a"] = df_scores.weight_a * df_scores.prediction_a
#             df_scores["score_b"] = df_scores.weight_b * df_scores.prediction_b
#             df_scores["score"] = df_scores.score_a + df_scores.score_b
#             summed_score = sum(df_scores.score)/2
#             match = re.match(".*/([A-Za-z0-9- ]+)_([A-Za-z0-9- ]+)_all.csv.bz2", distance_file)
#             drug_a = match.group(1)
#             drug_b = match.group(2)
#             n_paths = len(df_scores)
#             if n_paths != 0:
#                 mean_score = summed_score/n_paths
#             else:
#                 mean_score = summed_score
#             summed_scores.append(
#                 {
#                     "drug_a": drug_a,
#                     "drug_b": drug_b,
#                     "summed_score": summed_score,
#                     "mean_score": mean_score,
#                     "n_paths": n_paths
#                 }
#             )
# 
#         pd.DataFrame(summed_scores).to_csv(output.summed_scores, sep = "\t", index=False)
# 
# 
# rule calculate_drug_drug_distance:
#     input:
#         names = "work/chembl/pairs/name_conversion.csv",
#         from_to = "work/chembl/gene_interactions/{druga}_{drugb}.csv.bz2",
#         distance = "work/chembl/pairs/distance_pairs_possible_ends.csv"
#     output:
#         all_distances = "work/chembl/distances/{druga}_{drugb}_all.csv.bz2"
#     run:
#         name_dict = dict()
#         with open(input.names,"r") as f:
#             for line in f:
#                 gene_name, stringid, _ = line.split("\t")
#                 name_dict.update({gene_name: stringid})
# 
#         t = lambda x: name_dict[x]
#         df_from_to = pd.read_csv(input.from_to, sep="\t")
#         df_from_to = df_from_to[(df_from_to["gene_a"] != "Non-human") & (df_from_to["gene_b"] != "Non-human")]
#         df_from_to["gene_a"] = df_from_to["gene_a"].apply(t)
#         df_from_to["gene_b"] = df_from_to["gene_b"].apply(t)
#         df_distance = pd.read_csv(input.distance, sep="\t")
#         df_distance.columns.values[0:2] = ["gene_a", "gene_b"]
#         df_all_first = df_distance.merge(df_from_to, how="right", on=["gene_a", "gene_b"])
#         df_all_flip = df_all_first[(df_all_first.depth != df_all_first.depth)]
#         df_all_flip = df_all_flip[[
#             "gene_a",
#             "gene_b",
#             "target_a",
#             "target_b",
#             "prediction_a",
#             "prediction_b"
#         ]]
#         df_all_flip.columns = [
#             "gene_b",
#             "gene_a",
#             "target_b",
#             "target_a",
#             "prediction_b",
#             "prediction_a"
#         ]
#         df_all_flip = df_distance.merge(df_all_flip,how="right",on=["gene_a", "gene_b"])
#         df_all_first = df_all_first[~(df_all_first.depth != df_all_first.depth)]
#         df_all = pd.concat([df_all_first, df_all_flip])
#         df_all.to_csv(output.all_distances, sep="\t", index = False)
# 

rule calculate_gene_distance:
    input:
        all_unique = "work/HPO/pairs/distance_pairs.csv",
        genes = lambda wildcards: get_all_genes_at_depth(wildcards)
    output:
        distance = "work/HPO/pairs/distance_pairs_possible_ends.csv"
    run:
        def _get_depth(gene, target, depth, suffix):
            df_path = pd.read_csv(
                f"work/HPO/gene_path/{gene}_at_{depth}_shortest.csv.bz2",
                sep="\t"
            )
            number_hits = sum(df_path.iloc[:, -1] == target)
            number_of_paths = len(df_path.index)
            return {
                f"hits{suffix}": number_hits,
                f"total_paths{suffix}": number_of_paths,
                f"weight{suffix}": number_hits/number_of_paths
            }

        df_unique = pd.read_csv(input.all_unique, sep = "\t", header=None)
        expected_distance = []
        tot_rows = len(df_unique)
        for i, gene_a, gene_b, depth in df_unique.itertuples():
            row = {
                "from": gene_a,
                "to": gene_b,
                "depth": depth
            }
            row.update(
                _get_depth(gene_a, gene_b, depth, "_a")
            )
            row.update(
                _get_depth(gene_b, gene_a, depth, "_b")
            )
            expected_distance.append(row)
            print(f"{i}/{tot_rows} rows: {round(i/tot_rows, 3)*100} % done", end="\r")

        df = pd.DataFrame(expected_distance)
        df.to_csv(output.distance, sep = "\t", index=False)


rule guarantee_shortest_path:
    input:
        gene = "work/HPO/gene_path/{geneb}_at_{depth}.csv"
    output:
        gene = "work/HPO/gene_path/{geneb}_at_{depth}_shortest.csv.bz2"
    run:
        possible_paths = pd.read_csv(input.gene, sep = "\t")
        visited = possible_paths.iloc[:, :-1].unstack().unique()
        shortest_paths = possible_paths[~possible_paths.iloc[:, -1].isin(visited)]
        shortest_paths.to_csv(output.gene, sep = "\t", index=False)
# 
checkpoint get_possible_paths_from_genes:
    input:
        ppi_network_file = "data/9606.protein.links.v11.5.txt",
        unique_depth = "work/HPO/pairs/unique_depth.csv",
    output:
        directory("work/HPO/gene_path")
    shell:
        """
        mkdir -p {output}
        python src/HPO/get_all_paths.py {input.unique_depth} {input.ppi_network_file} {output}
        """


rule get_all_unique_depth:
    input:
        ppi_network_file = "data/9606.protein.links.v11.5.txt",
        distance_pairs = "work/HPO/pairs/distance_pairs.csv",
        names = "work/HPO/pairs/name_conversion.csv"
    output:
        all_unique_depths="work/HPO/pairs/unique_depth.csv"
    run:
        gp = GetGeneDistanceAtN(input.ppi_network_file)
        gp.load_string_identifiers(input.names)
        gp.load_previous_distances(input.distance_pairs)
        gp.get_all_at_depth(output.all_unique_depths)


n_parts = 100

rule distance_pairs_collect:
    input:
        distance_pairs = expand("work/HPO/pairs/distance_pairs_{i}.csv", i=range(n_parts))
    output:
        "work/HPO/pairs/distance_pairs.csv"
    shell:
        """
        awk 'FNR > 1' {input.distance_pairs} >> {output}
        """


rule get_gene_depth_pairs_HPO_part:
    """
    parallelled
    """
    input:
        ppi_network_file = "data/9606.protein.links.v11.5.txt",
        unique_interactions = "work/HPO/gene_interactions/unique_{i}.csv",
        names= "work/HPO/pairs/name_conversion.csv",
        missing= "work/HPO/pairs/missing.csv"
    output:
        distance_pairs = "work/HPO/pairs/distance_pairs_{i}.csv"

    run:
        gp = GetGeneDistanceAtN(input.ppi_network_file)
        gp.set_tmpdir(f"work/HPO/pairs/tmp_{wildcards.i}")
        gp.set_missing(input.missing)
        gp.load_combinations(input.unique_interactions)
        gp.load_string_identifiers(input.names)
        gp.calculated_distance_pairs()
        gp.write_gene_distances(output.distance_pairs)


rule get_unq_parts:
    input:
        unique = "work/HPO/gene_interactions/unique.csv"
    output:
        parts = expand("work/HPO/gene_interactions/unique_{i}.csv", i = range(n_parts))
    run:
        df_all = pd.read_csv(input.unique, sep = "\t")
        n_parts = len(output.parts)
        nrows = len(df_all.index)
        step = int(nrows/n_parts)
        part_row = range(0, nrows, step)
        for i, current in enumerate(output.parts):
            df_all.iloc[range(part_row[i], part_row[i+1])].to_csv(current, index = False, sep="\t")

rule get_name_df_hpo:
    input:
        ppi_network_file = "data/9606.protein.links.v11.5.txt",
        unique ="work/HPO/gene_interactions/unique.csv"
    output:
        names = "work/HPO/pairs/name_conversion.csv",
        missing = "work/HPO/pairs/missing.csv"
    run:

        gp = GetGeneDistanceAtN(input.ppi_network_file)
        gp.set_missing(output.missing)
        gp.get_format_string_identifiers(input.unique)
        gp.write_name_df(output.names)

rule get_unique_gene_interaction_hpo:
    input:
        all_dfs = [f"work/HPO/gene_interactions/{hpo_a}_{hpo_b}.csv.bz2"
                   for hpo_a, hpo_b in hpo_combinations]
    output:
        all_unique = "work/HPO/gene_interactions/unique.csv"
    run:

        unique_set = set()
        n = len(input.all_dfs)
        for interaction_file_name  in input.all_dfs:
            with bz2.open(interaction_file_name, "r") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.decode("utf-8").strip().split("\t")
                    pair = tuple(sorted([line[0], line[1]]))
                    if pair not in unique_set:
                        unique_set.update({pair})

        with open(output.all_unique, "w") as w:
            w.write("gene_a\tgene_b\n")
            for a, b in unique_set:
                w.write(f"{a}\t{b}\n")




rule get_drug_drug_gene_interaction_hpo:
    input:
        hpo_a = "work/HPO/interaction/{hpo_a}.csv",
        hpo_b= "work/HPO/interaction/{hpo_b}.csv",
    output:
        gene_interactions = "work/HPO/gene_interactions/{hpo_a}_{hpo_b}.csv.bz2"
    run:
        df_a = pd.read_csv(input.hpo_a, sep="\t")
        df_b = pd.read_csv(input.hpo_b, sep="\t")
        df_combinations = df_a.merge(df_b, how="cross", suffixes=("_a", "_b"))
        df_combinations.to_csv(output.gene_interactions, sep = "\t", index=False)

rule interacting_genes_hpo:
    input:
        hpo_csv = "data/MatchedHPO.csv"
    output:
        expand("work/HPO/interaction/{hpo}.csv", hpo = all_hpos)
    run:
        ghpo = GetHPOGenes(input.hpo_csv, "work/HPO/interaction")
        ghpo.query_hpo()


