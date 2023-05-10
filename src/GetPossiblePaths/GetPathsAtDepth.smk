import pandas as pd
import bz2

checkpoint get_possible_paths_from_genes:
    params:
        max_depth = config["max_depth"]
    input:
        ppi_network_file = "data/9606.protein.links.v11.5.txt",
        unique_genes = "work/{project}/gene_interactions/unique_genes.csv"
    output:
        directory("work/{project}/gene_path")
    singularity:
        "data/igraph.simg"
    shell:
        """
        mkdir -p {output}
        python src/GetPossiblePaths/GetPossiblePaths.py {input.unique_genes} {input.ppi_network_file} {output} {params.max_depth}
        """


rule guarantee_shortest_path:
    input:
        gene = "work/{project}/gene_path/{gene}_at_{depth}.csv"
    output:
        gene = "work/{project}/gene_path/{gene}_at_{depth}_shortest.csv.bz2"
    run:
        possible_paths = pd.read_csv(input.gene, sep = "\t")
        visited = possible_paths.iloc[:, :-1].unstack().unique()
        shortest_paths = possible_paths[~possible_paths.iloc[:, -1].isin(visited)]
        shortest_paths.to_csv(output.gene, sep = "\t", index=False)


rule possible_endpoints:
    input:
        full_paths = "work/{project}/gene_path/{gene}_at_{depth}_shortest.csv.bz2"
    output:
        summed_paths = "work/{project}/gene_path/{gene}_at_{depth}_sum.csv.bz2"
    run:
        gene_count = dict()
        with bz2.open(input.full_paths, "r") as f:
            header = True
            for line in f:
                if header:
                    header = False
                else:
                    line = line.decode("utf-8").strip().split()
                    path_end = line[-1]
                    if path_end in gene_count:
                        gene_count[path_end] += 1
                    else:
                        gene_count[path_end] = 1

        with bz2.open(output.summed_paths, "wt") as w:
            w.write("Gene\tCount\n")
            for gene, count in gene_count.items():
                w.write(f"{gene}\t{count}\n")

