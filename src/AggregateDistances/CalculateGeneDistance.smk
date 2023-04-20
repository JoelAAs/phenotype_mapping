import pandas as pd
import bz2
import os


#### Input Functions
def get_all_genes_at_depth(wc):

    all_needed_files_exists = all(
        map(
            os.path.exists,
        [
            f"work/{wc.project}/gene_path",
            f"work/{wc.project}/depth/unique.csv"]))
    if not all_needed_files_exists:
        re_eval = checkpoints.get_possible_paths_from_genes.get(**wc)
        return re_eval
    else:
        df_distance_pairs = pd.read_csv(f"work/{wc.project}/depth/unique.csv", sep = "\t")
        expected_distance = []
        for i, gene, depth in df_distance_pairs.itertuples():
            print(gene)
            expected_distance.append(
                f"work/{wc.project}/gene_path/{gene}_at_{depth}_sum.csv.bz2"
            )

        return expected_distance


#### Rules
rule get_list_of_all_gene_depths:
    input:
        genes = lambda wildcards: get_all_genes_at_depth(wildcards)
    output:
        expected_output = "work/{project}/checkpoint/all_genes_at_depth.txt"
    run:
        with open(output.expected_output, "w") as w:
            for gene in input.genes:
                w.write(f"{gene}\n")


rule calculate_gene_distance:
    input:
        all_unique = "work/{project}/depth/all_depths.csv",
        gene_files = "work/{project}/checkpoint/all_genes_at_depth.txt"
    output:
        distance = "work/{project}/gene_distance/all_unique_gene_distance.csv"
    run:
        def _get_depth(gene, target, depth):
            with bz2.open(f"work/{wildcards.project}/gene_path/{gene}_at_{depth}_sum.csv.bz2", "r") as f:
                number_of_paths = 0
                number_hits = 0
                all_lines = f.readlines()
                for gene_path_line in all_lines[1:]:
                    gene_path_line = gene_path_line.decode("utf-8").strip().split()
                    if gene_path_line[0] == target:
                        number_hits = int(gene_path_line[1])
                    number_of_paths += int(gene_path_line[1])

            return [
                number_hits,
                number_of_paths,
                number_hits/number_of_paths
                ]

        df_unique = pd.read_csv(input.all_unique, sep = "\t")
        tot_rows = len(df_unique)
        with open(output.distance, "w") as w:
            w.write("gene_a\tgene_b\tdepth\thit_a\ttot_path_a\tweight_a\thit_b\ttot_path_b\tweight_b\n")
            for i, gene_a, gene_b, depth in df_unique.itertuples():
                hit_b, tot_path_b, weight_b = _get_depth(gene_b, gene_a, depth)
                hit_a, tot_path_a, weight_a = _get_depth(gene_a, gene_b, depth)
                print(f"{i}/{tot_rows} rows: {round(i/tot_rows, 3)*100} % done", end="\r")
                w.write(f"{gene_a}\t{gene_b}\t{depth}\t{hit_a}\t{tot_path_a}\t{weight_a}\t{hit_b}\t{tot_path_b}\t{weight_b}\n")


