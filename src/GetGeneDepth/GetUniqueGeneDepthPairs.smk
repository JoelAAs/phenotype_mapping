

rule get_unique_gene_depth_pairs:
    input:
        all_depths = "work/{project}/depth/all_depths.csv"
    output:
        unique_gene_depths = "work/{project}/depth/unique.csv"
    run:
        unique_depths = set()
        header = True
        with open(input.all_depths, "r") as f:
            for line in f:
                if header:
                    header = False
                else:
                    gene_a, gene_b, depth = line.strip().split("\t")
                    unique_depths.update({(gene_a, depth), (gene_b, depth)})

        with open(output.unique_gene_depths, "w") as w:
            w.write(f"gene\tdepth\n")
            for gene, depth in unique_depths:
                w.write(f"{gene}\t{depth}\n")

