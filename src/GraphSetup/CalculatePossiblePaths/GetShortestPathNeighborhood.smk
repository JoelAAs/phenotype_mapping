import pandas as pd
import bz2
import os
import pathlib
import glob

rule unique_genes:
    output:
        unique_genes = expand("work/{{project}}/unique_genes_{n}.csv", n = range(config["n_path_batches"]))
    run:
        for n in range(config["n_path_batches"]):
            outfile = f"work/{wildcards.project}/unique_genes_{n}.csv"
            with open(outfile, "w") as w:
                w.write("Gene\n")
                for gene in config["path_batches"][n]:
                    w.write(f"{gene}\n")


rule get_possible_paths_from_genes:
    params:
        max_depth = config["max_depth"],
        ppi_network_file= config["ppi_file"]
    input:
        unique_genes = "work/{project}/unique_genes_{n}.csv"
    output:
        folder =  directory("work/{project}/paths_{n}")
    singularity:
        "data/igraph.simg"
    shell:
        """
        mkdir -p {output.folder}
        python3 src/GraphSetup/CalculatePossiblePaths/possible_paths_from_a.py \
            {input.unique_genes} \
            {params.ppi_network_file} \
            {output.folder} \
            {params.max_depth}
        """


rule check_for_all_paths:
    """
    Delete batch_folders files if changes are made to code in possible_paths_from_a.py
    """
    input:
        folders = expand("work/{{project}}/paths_{n}", n = range(config["n_path_batches"]))
    output:
        all_genes = expand(
            "work/{{project}}/paths/{gene}_at_{depth}.csv",
            gene=config["all_unique_genes"],
            depth=range(1,config["max_depth"] + 1)
        )
    run:
        for batch_folder in input.folders:
            files = glob.glob(batch_folder + "/*")
            for path_file in files:
                path_filename = os.path.basename(path_file)
                abs_path_file = os.path.abspath(path_file)
                ln_file = f"work/{wildcards.project}/paths/{path_filename}"
                shell(f"ln -s {abs_path_file} {ln_file}")


rule guarantee_shortest_path:
    input:
        gene = "work/{project}/paths/{gene}_at_{depth}.csv"
    output:
        gene = "work/{project}/paths/{gene}_at_{depth}_shortest.csv.bz2"
    run:
        possible_paths = pd.read_csv(input.gene, sep = "\t")
        visited = possible_paths.iloc[:, :-1].unstack().unique()
        shortest_paths = possible_paths[~possible_paths.iloc[:, -1].isin(visited)]
        shortest_paths.to_csv(output.gene, sep = "\t", index=False)


rule genes_in_path_neighborhood:
    """
    Calculates probabilities of picking gene out of all shortest paths originating from gene
    """
    input:
        depth_genes = expand("work/{{project}}/paths/{{gene}}_at_{depth}_shortest.csv.bz2",
            depth = range(1, config["max_depth"] + 1))
    output:
        genes_in_paths = "work/{project}/neighborhood/{gene}_p_gene.csv.bz2"
    run:
        gene_count = dict()
        n_paths = 0
        for at_depth in input.depth_genes:
            with bz2.open(at_depth, "r") as f:
                header = True
                for line in f:
                    if header:
                        header = False
                    else:
                        n_paths += 1
                        line = line.decode("utf-8").strip().split()
                        for gene in line:
                            if gene in gene_count:
                                gene_count[gene] += 1/len(line)  ## sum(gene) over P(m|gene, b)p(gene|path)
                            else:
                                gene_count[gene] = 1/len(line)

        with bz2.open(output.genes_in_paths, "wt") as w:
            w.write("Gene\tP_gene\n")
            for gene, p_gene_path in gene_count.items():
                w.write(f"{gene}\t{p_gene_path/n_paths}\n") ## P(m) = P(m|path)P(path)
