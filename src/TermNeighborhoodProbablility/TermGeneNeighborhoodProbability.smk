import bz2
import pandas as pd
## input function

def get_genes_in_query(wc):
    genes = config[wc["term2"]]
    expected = [
        f"work/{wc.project}/term_gene_distance/{wc.term2}_{gene}_distance.csv" for gene in genes
        ]
    return expected

## rules
rule term_gene_probability:
    input:
        term_geneset = "input/{project}/{term}.csv",
        gene_neighborhood = "work/{project}/neighborhood/{gene}_p_gene.csv.bz2"
    output:
        term_prob = "work/{project}/term_gene_distance/{term}_{gene}_distance.csv"
    run:
        probability_dict = {}
        with bz2.open(input.gene_neighborhood, "r") as f:
            header = True
            for line in f:
                if header:
                    header = False
                else:
                    gene, prob = line.decode("utf-8").strip().split("\t")
                    probability_dict[gene] = prob

        with open(output.term_prob, "w") as w:
            w.write("gene\tprobability\n")
            with open(input.term_geneset, "r") as f:
                term_geneset = [l.strip() for l in f]
                for term_gene in term_geneset:
                    if term_gene in probability_dict:
                        w.write(f"{term_gene}\t{probability_dict[term_gene]}\n")
                    else:
                        w.write(f"{term_gene}\t{0}\n")

rule full_term_to_term_matrix:
    input:
        neighborhood_term_geneset = "input/{project}/{term1}.csv",
        query_term_geneset = "input/{project}/{term2}.csv",
        term_distances = get_genes_in_query
    output:
        full_matrix = "work/{project}/term_gene_distance/full_matrix_{term2}_in_{term1}.csv"
    run:
        df_list = []
        first = True
        for term_gene_probability_file in input.term_distances:
            gene_neighborhood = term_gene_probability_file.split("-")[-2]
            tgpf_df = pd.read_csv(term_gene_probability_file, sep = "\t")
            tgpf_df = tgpf_df.rename({"gene": "gene", "probability": gene_neighborhood})
            if not first:
                full_df = full_df.join(tgpf_df, on="gene")
            else:
                full_df = tgpf_df
                first = False
        full_df.to_csv(output.full_matrix, sep = "\t")

