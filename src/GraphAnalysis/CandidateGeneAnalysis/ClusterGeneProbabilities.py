import pandas as pd
import glob
import networkx as nx


def get_occurrence_of_genes_in_terms(df_cluster_definition_file, cluster_notation, projects):
    clusters_df = pd.read_csv(df_cluster_definition_file, sep="\t")
    clusters_df = clusters_df[clusters_df["Cluster"] == int(cluster_notation)]

    all_genes = []
    n_terms = len(clusters_df["Node"])
    for i, row in clusters_df.iterrows():
        node = row["Node"]
        for project in projects:
            geneset_file = glob.glob(f"input/{project}/{node}.csv")
            if len(geneset_file) > 1:
                raise IOError(f"Multiple hits for inputfile for {node} in {project}")
            if geneset_file:
                with open(geneset_file[0], "r") as f:
                    genes = [l.strip() for l in f][1:]
                    all_genes += genes

    all_genes = pd.Series(all_genes).value_counts()
    count_df = pd.DataFrame(
        {"gene": all_genes.index,
         "count": all_genes.values}
    )

    return count_df, n_terms


def get_centrality_of_cluster_genes(ppi_file, probabilities_df, top_fraction=0.6):
    edge_list_df = pd.read_csv(ppi_file, sep=" ")
    edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
    edge_list_df = edge_list_df[edge_list_df["combined_score"] > top_fraction]
    G = nx.from_pandas_edgelist(
        edge_list_df,
        "protein1",
        "protein2",
        create_using=nx.Graph
    )
    centrality_dict = nx.degree_centrality(G)
    probability_df = pd.read_csv(probabilities_df, sep="\t")
    probability_df["centrality"] = probability_df.apply(
        lambda row: centrality_dict[row["gene"]], axis=1)
    return probability_df
