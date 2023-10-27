import pandas as pd
import bz2


rule get_possible_paths_from_genes:
    input:
        ppi_network_file = "data/9606.protein.links.v11.5.txt",
        unique_genes = "work/{project}/gene_interactions/unique_genes.csv"
    output:
        all_genes = expand(
            "work/{project}/gene_path/{gene}_at_{depth}.csv",
            project = config["project_name"],
            gene = config["all_unique_genes"],
            depth = range(config["max_depth"]) ## TODO length
        )
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


rule genes_in_path_neighborhood:
    input:
        depth_genes = expand("work/{project}/gene_path/{gene}_at_{depth}_shortest.csv.bz2",
            depth = range(1, config["max_depth"]))
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
                                gene_count[gene] += 1/len(line)  ## P(gene|path)
                            else:
                                gene_count[gene] = 1/len(line)

        with bz2.open(output.genes_in_paths, "wt") as w:
            w.write("Gene\tP_gene\n")
            for gene, p_gene_path in gene_count.items():
                w.write(f"{gene}\t{p_gene_path/n_paths}\n") ## P(gene) = P(gene|path)P(path)
