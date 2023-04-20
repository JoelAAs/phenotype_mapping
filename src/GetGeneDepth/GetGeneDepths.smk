from CalculateGeneDistance import CalculateGeneDistance
import pandas as pd


wildcard_constraints:
    part="[0-9]+"

rule get_gene_pair_parts:
    input:
        unique_pairs = "work/{project}/gene_interactions/unique.csv"
    output:
        parts = expand("work/{{project}}/gene_interactions/unique_{part}.csv", part = range(config["depth_parts"]))
    run:
        df_all = pd.read_csv(input.unique_pairs,sep="\t")
        n_parts = len(output.parts)
        nrows = len(df_all.index)
        step = int(nrows/n_parts)
        part_row = range(0, nrows, step-1)
        for i, current in enumerate(output.parts):
            print(part_row)
            print(df_all)
            df_all.iloc[range(part_row[i], part_row[i + 1])].to_csv(current,index=False,sep="\t")

rule get_gene_pair_depths:
    params:
        ppi_location = config["ppi_file"],
    input:
        unique_pairs = "work/{project}/gene_interactions/unique_{part}.csv"
    output:
        depth_part = "work/{project}/depth/{part}.csv"
    run:
        cdg = CalculateGeneDistance(params.ppi_location)
        cdg.calculate_depth_of_pairs(input.unique_pairs, output.depth_part)


rule concat_depth_parts:
    input:
        expand("work/{{project}}/depth/{part}.csv", part=range(config["depth_parts"]))
    output:
        "work/{project}/depth/all_depths.csv"
    shell:
        """
        echo "gene_a\tgene_b\tdistance" > {output}
        awk 'FNR > 1' {input} >> {output}
        """
