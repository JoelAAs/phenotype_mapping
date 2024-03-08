import math
import os
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import lognorm, gumbel_r, weibull_min
def get_components(wildcards):
    comp_check = checkpoints.get_connected_components.get(**wildcards).output[0]
    comps, = glob_wildcards(os.path.join(comp_check, "Component_{c}.csv"))

    expected = f"work/full-drugbank-benchmark/group-quant/{wildcards.group}_connected_components_enrichment/Component_{{component}}_KEGG.png"
    return expand(expected, component = comps)

rule get_CDF_of_genes:
    params:
        ppi_file = "data/9606.protein.links_above_700.v11.5.txt",
        limit = 0.7,
        n_permuts = config["n_permutation_combinations"]
    input:
        gene_probabilities = "data/full_drugbank_gene_probabilites/{group}_gene_probabilites.csv",
        set_permutations = "work/full-drugbank-benchmark/candidate_genes/group_gene_permut/{group}_permut.csv"
    output:
        quant = "work/full-drugbank-benchmark/group-quant/{group}_quant.csv",

    run:
        set_permut_dict = dict()
        with open(input.set_permutations, "r") as w:
            lines = w.readlines()[1:]
            for l in lines:
                gene, prob, i = l.strip().split("\t")
                prob = float(prob)
                if gene not in set_permut_dict:
                    set_permut_dict[gene] = np.full(params.n_permuts, fill_value=0.0)
                set_permut_dict[gene][int(i)] = prob

        all_data = np.full(len(set_permut_dict)*params.n_permuts, fill_value=0.0)
        for i, (gene, probs) in enumerate(set_permut_dict.items()):
            all_data[i*params.n_permuts:(i+1)*params.n_permuts] = probs

        mean_loc, mean_shape = gumbel_r.fit(all_data)

        # PPI Degree conf
        edge_list_df = pd.read_csv(params.ppi_file,sep=" ")
        edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
        edge_list_df = edge_list_df[edge_list_df["combined_score"] >= params.limit]
        G = nx.from_pandas_edgelist(
            edge_list_df,
            "protein1",
            "protein2"
        )

        degree_dict = {}
        for node, c_degree in G.degree():
            if c_degree in degree_dict:
                degree_dict[c_degree].append(node)
            else:
                degree_dict[c_degree] = [node]
        node_degree_dict = {node: degree for node, degree in G.degree()}

        mean_var_dict = dict()
        def get_sample_mean_var(gene):
            if gene not in set_permut_dict:
                return mean_loc, mean_shape, 1
            else:
                probs = set_permut_dict[gene]
                s_loc, s_scale = gumbel_r.fit(probs)
            return s_loc, s_scale, 0

        with open(input.gene_probabilities, "r") as f:
            with open(output.quant, "w") as w:
                w.write("Gene\ty_probability\tquant\tdegree\tnot_reached\ts_loc\ts_scale\n")
                lines = f.readlines()[1:]
                genes = np.full(len(lines), fill_value="",dtype=object)
                for i, l in enumerate(lines):
                    gene, prob = l.strip().split("\t")
                    genes[i] = gene
                    prob = float(prob)
                    if gene in node_degree_dict:
                        s_loc, s_scale, missing = get_sample_mean_var(gene)
                        quant = gumbel_r.cdf(prob, loc=s_loc, scale=s_scale)
                        w.write(f"{gene}\t{prob}\t{quant}\t{node_degree_dict[gene]}\t{missing}\t{s_loc}\t{s_scale}\n")

                not_reached = [
                    gene_node for gene_node in node_degree_dict.keys() if gene_node not in genes
                ]
                for gene_left in not_reached:
                    s_loc, s_scale, missing = get_sample_mean_var(gene_left)
                    quant = gumbel_r.cdf(0, loc=s_loc, scale=s_scale)
                    w.write(f"{gene_left}\t{0}\t{quant}\t{node_degree_dict[gene_left]}\t{missing}\t{s_loc}\t{s_scale}\n")


rule get_top_values:
    params:
        n_top = 500,
        string_names = "data/stringdb/9606.protein.info.v11.5.txt",
        keep_input = True
    input:
        z_scores = "work/full-drugbank-benchmark/group-quant/{group}_quant.csv",
        unique_input = "data/full_drugbank_gene_probabilites/{group}_unique_genes.csv"
    output:
        top_z_scores = "work/full-drugbank-benchmark/group-quant/{group}_top.csv"
    run:
        z_scores = pd.read_csv(input.z_scores, sep="\t")
        with open(input.unique_input, "r") as f:
            input_genes = [g.strip() for g in f.readlines()[1:]]

        string_df = pd.read_csv(params.string_names, sep = "\t")[["string_protein_id", "preferred_name"]]
        z_scores = z_scores.merge(string_df, left_on="Gene", right_on="string_protein_id", how="left")
        if not params.keep_input:
            z_scores = z_scores[~z_scores["Gene"].isin(input_genes)]
        del z_scores["string_protein_id"]
        z_scores_top = z_scores.sort_values(by="quant", ascending=False)[:(params.n_top + len(input_genes))]
        z_scores_top.to_csv(output.top_z_scores, index=False, sep = "\t")


checkpoint get_connected_components:
    params:
        min_subgraph_size = 3,
        ppi_file = "data/9606.protein.links.v11.5.txt",
        limit = 0.7
    input:
        top_z_scores = "work/full-drugbank-benchmark/group-quant/{group}_top.csv",
    output:
        output_dir = directory("work/full-drugbank-benchmark/group-quant/{group}_top_connected_components")
    run:

        edge_list_df = pd.read_csv(params.ppi_file,sep=" ")
        edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
        edge_list_df = edge_list_df[edge_list_df["combined_score"] >= params.limit]
        G = nx.from_pandas_edgelist(
            edge_list_df,
            "protein1",
            "protein2"
        )
        top_df = pd.read_csv(input.top_z_scores, sep="\t")
        keep_genes = set(top_df["Gene"])
        remove_nodes = [n for n in G.nodes() if n not in keep_genes]
        G.remove_nodes_from(remove_nodes)

        subgraph_n = 0
        os.mkdir(output.output_dir)
        outfile = output.output_dir + "/Component_{component}.csv"
        for cc in sorted(nx.connected_components(G), key=len, reverse=True):
            if len(cc) > params.min_subgraph_size:
                with open(outfile.format(component=subgraph_n), "w") as w:
                    subgraph_n += 1
                    w.write("Gene\n")
                    for gene in cc:
                        w.write(gene + "\n")

rule remove_input_and_translate_to_entrez:
    input:
        unique_input = "data/full_drugbank_gene_probabilites/{group}_unique_genes.csv",
        entrez= "data/ncbi/entrez.csv",
        string_id="data/stringdb/9606.protein.info.v11.5.txt",
        components = "work/full-drugbank-benchmark/group-quant/{group}_top_connected_components/Component_{component}.csv"
    output:
        translated = "work/full-drugbank-benchmark/group-quant/{group}_connected_components_enrichment/Component_{component}.csv",
    run:
        with open(input.unique_input, "r") as f:
            input_genes = [l.strip() for l in f.readlines()[1:]]

        components_genes = pd.read_csv(input.components, sep = "\t")
        entrez_df = pd.read_csv(input.entrez, sep="\t")
        string_df = pd.read_csv(input.string_id, sep = "\t")

        string_genes = components_genes.merge(
            string_df,
            left_on = "Gene",
            right_on="string_protein_id",
            how="left")

        entrez_genes = string_genes.merge(
            entrez_df,
            left_on="preferred_name",
            right_on="gene_name",
            how="left")

        entrez_genes[["gene_name", "string_protein_id", "entrez"]].to_csv(output.translated, sep = "\t", index=False)

rule component_KEGG_enrich:
    input:
        translated = "work/full-drugbank-benchmark/group-quant/{group}_connected_components_enrichment/Component_{component}.csv"
    output:
        results = "work/full-drugbank-benchmark/group-quant/{group}_connected_components_enrichment/Component_{component}_KEGG.csv",
    shell:
        """
        Rscript src/GraphAnalysis/CandidateGeneAnalysis/Enrichment.R {input} {output}
        """


rule enrichment_done_components:
    input:
        get_components
    output:
        done = "work/full-drugbank-benchmark/group-quant/{group}_connected_components_enrichment/done.csv"
    shell:
        """
        touch {output.done}
        """


rule get_all_shortest_distances:
    input:
        quants = "work/full-drugbank-benchmark/group-quant/{group}_quant.csv",
        unique_input = "data/full_drugbank_gene_probabilites/{group}_unique_genes.csv",
        ppi_file = "data/9606.protein.links_above_700.v11.5.txt"
    output:
        distances = "work/full-drugbank-benchmark/group-quant/shortest_distances/{group}.csv"
    run:
        edge_list_df = pd.read_csv(input.ppi_file,sep=" ")
        edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
        edge_list_df = edge_list_df[edge_list_df["combined_score"] >= 0.7]
        G = nx.from_pandas_edgelist(
            edge_list_df,
            "protein1",
            "protein2"
        )
        input_genes_df = pd.read_csv(input.unique_input,sep="\t")
        genes_input = input_genes_df["gene"].tolist()
        prob_genes_df = pd.read_csv(input.quants, sep = "\t")
        genes_to = prob_genes_df["Gene"].tolist()
        with open(output.distances, "w") as w:
            w.write("input_gene\ttarget_gene\tdistance\n")
            for i, gene_input in enumerate(genes_input):
                print(f"{i} of {len(genes_input)}")
                for gene_to in genes_to:
                    try:
                        d = nx.shortest_path_length(
                            G,
                            source=gene_input,
                            target=gene_to
                        )
                        w.write(f"{gene_input}\t{gene_to}\t{d}\n")
                    except nx.NetworkXNoPath as e:
                        pass

rule get_average_shortest_distances:
    input:
        distances = "work/full-drugbank-benchmark/group-quant/shortest_distances/{group}.csv",
        counts = "data/full_drugbank_gene_probabilites/{group}_count_unique_genes.csv"
    output:
        mean_dist = "work/full-drugbank-benchmark/group-quant/shortest_distances/{group}_proportinal.csv"
    run:
        mean_dist_df = pd.read_csv(input.counts, sep="\t")
        mean_dist_df["proportion"] = mean_dist_df.Count/mean_dist_df.Count.sum()
        distances_df = pd.read_csv(input.distances, sep="\t")
        distances_count_df = distances_df.merge(mean_dist_df, left_on="input_gene", right_on="Gene")
        distances_count_df["proportional_distance"] = distances_count_df["proportion"]*distances_count_df["distance"]
        proportional_distances = distances_count_df.groupby("target_gene", as_index=False)["proportional_distance"].sum()
        proportional_distances.to_csv(output.mean_dist, sep="\t", index=False)

