import pandas as pd
from Chembl.GetDistancePairs import GetGeneDistanceAtN
from itertools import product
import re
import os

def  get_all_medication(file_csv):
    df = pd.read_csv(file_csv, sep = "\t")
    return list(df.ENG)

drugs = get_all_medication("data/Drugs.csv")

wildcard_constraints:
    drug="[A-Za-z1-9. \-()]+"

### Input functions

def possible_gene_combinations(wildcards):
    re_eval = checkpoints.get_possible_paths_from_genes.get()
    name_dict = dict()
    with open("work/chembl/pairs/name_conversion.csv", "r") as f:
        for line in f:
            gene_name, stringid, _ = line.split("\t")
            name_dict.update({gene_name: stringid})

    depth_dict = dict()
    with open("work/chembl/pairs/distance_pairs.csv","r") as f:
        for line in f:
            genea, geneb, depth = line.strip().split("\t")
            dict_key = "".join(sorted([genea, geneb]))
            depth_dict.update({dict_key: depth})

    df_a = pd.read_csv(
        f"work/chembl/predicted_interaction/annotated/{wildcards.druga}.csv",
        sep = "\t"
    )
    df_a = df_a[~df_a.gene == "Non-human"]
    genes_a = [name_dict[gene] for gene in list(df_a["target"])]
    df_b = pd.read_csv(
        f"work/chembl/predicted_interaction/annotated/{wildcards.drugb}.csv",
        sep = "\t"
    )
    df_b = df_b[~df_b.gene == "Non-human"]
    genes_b = [name_dict[gene] for gene in list(df_b["target"])]
    gene_combinations = list(product(genes_a, genes_b))

    expected_input = [""] *len(gene_combinations)
    for i, a, b in enumerate(gene_combinations):
        expected_input[i] = "work/chembl/gene_path/{genea}_{geneb}_distance_at_{depth}.csv".format(
            genea=a, geneb=b, depth=depth_dict["".join(sorted([a,b]))]
        )

    return expected_input


rule calculate_drug_drug_distance:
    input:
        names = "work/chembl/pairs/name_conversion.csv",
        shortest_path ="work/chembl/pairs/distance_pairs.csv",
        drug_a = "work/chembl/predicted_interaction/annotated/{druga}.csv",
        drug_b = "work/chembl/predicted_interaction/annotated/{drugb}.csv",
        score_files = lambda wildcards: possible_gene_combinations(wildcards)
    output:
        all_distances = "work/chembl/distances/{druga}_{drugb}_all.csv"
    run:
        name_dict = dict()
        with open(input.names,"r") as f:
            for line in f:
                gene_name, stringid, _ = line.split("\t")
                name_dict.update({gene_name: stringid})
        drug_a = dict()
        with open(input.drug_a, "r") as f:
            lines = f.readlines()[1:]
            for line in lines:
                _, gene, weight = line.strip().split("\t")
                drug_a.update({name_dict[gene]: float(weight)})

        drug_b = dict()
        with open(input.drug_b, "r") as f:
            lines = f.readlines()[1:]
            for line in lines:
                _, gene, weight = line.strip().split("\t")
                drug_b.update({name_dict[gene]: float(weight)})

        with open(output.all_distances, "w") as w:
            w.write("from\tto\tdistance\tweight\tscore\n")
            for combination_score_file in input.score_files:
                matching_group = "([A-Za-z0-9\\.]+)"
                matches = re.search(
                    f"work/chembl/gene_path/{matching_group}_{matching_group}_distance_at_\d+.csv",
                    combination_score_file)
                gene_a =  matches.group(1)
                gene_b = matches.group(1)
                with(combination_score_file, "r") as f:
                    lines = f.readlines()
                    a, b, distancea2b = lines[1].strip().split("\t")
                    _, _, distanceb2a = lines[2].strip().split("\t")
                    w.write(lines[1].strip() + f"\t{drug_a[gene_a]}\t{drug_a[gene_a]*float(distancea2b)}\n")
                    w.write(lines[2].strip() + f"\t{drug_b[gene_b]}\t{drug_b[gene_b]*float(distanceb2a)}\n")


rule calculate_gene_distance:
    input:
        gene_a = "work/chembl/gene_path/{genea}_at_{depth}_shortest.csv",
        gene_b = "work/chembl/gene_path/{geneb}_at_{depth}_shortest.csv"
    output:
        distance = "work/chembl/gene_path/{genea}_{geneb}_distance_at_{depth}.csv"
    run:
        if gene_a == gene_b:
            with open(output.distance,"w") as f:
                f.write("From\tTo\tdistance\n")
                f.write(f"{wildcards.genea}\t{wildcards.geneb}\t1\n")
                f.write(f"{wildcards.geneb}\t{wildcards.genea}\t1\n")
        else:
            genea_df = pd.read_csv(input.gene_a)
            geneb_df = pd.read_csv(input.gene_b)

            a2b_name_check = genea_df.iloc[:, -1:].eq(wildcards.geneb)
            a2b = a2b_name_check.sum()[0]/len(a2b_name_check)

            b2a_name_check = genea_df.iloc[:, -1:].eq(wildcards.genea)
            b2a = b2a_name_check.sum()[0]/len(b2a_name_check)
            with open(output.distance,"w") as f:
                f.write("From\tTo\tdistance\n")
                f.write(f"{wildcards.genea}\t{wildcards.geneb}\t{a2b}\n")
                f.write(f"{wildcards.geneb}\t{wildcards.genea}\t{b2a}\n")


rule guarantee_shortest_path:
    input:
        gene = "work/chembl/gene_path/{geneb}_at_{depth}.csv"
    output:
        gene = "work/chembl/gene_path/{geneb}_at_{depth}_shortest.csv"
    run:
        possible_paths = pd.read_csv(input.gene, sep = "\t")
        visited = possible_paths.iloc[:, :-1].unstack().unique()
        shortest_paths = possible_paths[~possible_paths.iloc[:, -1].isin(visited)]
        shortest_paths.to_csv(output.gene, sep = "\t", index=False)

checkpoint get_possible_paths_from_genes:
    input:
        ppi_network_file = "data/9606.protein.links.v11.5.txt",
        distance_pairs = "work/chembl/pairs/distance_pairs.csv",
        names = "work/chembl/pairs/name_conversion.csv"
    output:
        directory("work/chembl/gene_path")
    run:
        gp = GetGeneDistanceAtN(input.ppi_network_file, output[0], "", [])
        gp.load_previous_distances(input.distance_pairs, input.names)
        gp.get_all_neighbours_at_depth()


rule get_gene_depth_pairs:
    """
    parallelled
    """
    input:
        ppi_network_file = "data/9606.protein.links.v11.5.txt",
        genes = expand("work/predicted_interaction/chembl/annotated/{drug}.csv", drug = drugs)
    output:
        distance_pairs = "work/chembl/pairs/distance_pairs.csv",
        names = "work/chembl/pairs/name_conversion.csv"
    threads:
        workflow.cores
    run:
        gp = GetGeneDistanceAtN(input.ppi_network_file, "", output.names, *input.genes)
        gp.calculated_distance_pairs()
        gp.write_gene_distances(output.distance_pairs)

rule get_genes:
    input:
        "work/chembl/predicted_interaction/{drug}.csv"
    output:
        "work/chembl/predicted_interaction/annotated/{drug}.csv"
    shell:
        """
        python src/Chembl/ChemblExtrainfo.py "{input}" "{output}"
        """

rule predict_interaction:
    output:
        "work/chembl/predicted_interaction/{drug}.csv"
    conda:
        "multitask-network"
    shell:
        """
        python src/Chembl/GetPredictedInteraction.py "{wildcards.drug}" "{output}"
        """



