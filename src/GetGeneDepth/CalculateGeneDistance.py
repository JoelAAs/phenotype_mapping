import pandas as pd
import igraph as ig


class CalculateGeneDistance:
    def __init__(self, ppi_network_file, min_score=0.6):
        self.graph = self.setup_network(ppi_network_file, min_score)

    @staticmethod
    def setup_network(ppi_network_file, min_score):
        edge_list_df = pd.read_csv(ppi_network_file, sep=" ")
        edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
        edge_list_df = edge_list_df[edge_list_df["combined_score"] > min_score]
        tuple_list = [tuple(x) for x in edge_list_df.values]
        ppi_graph = ig.Graph.TupleList(
            tuple_list,
            directed=False,
            weights=True)
        return ppi_graph

    def calculate_depth_of_pairs(self, combination_file, output_depths):
        with open(output_depths, "w") as w:
            w.write("gene_a\tgene_b\tdistance\n")
            with open(combination_file, "r") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    gene_a, gene_b = line.strip().split("\t")
                    depth = len(self.graph.get_all_shortest_paths(gene_a, gene_b)[0])
                    w.write(f"{gene_a}\t{gene_b}\t{depth}\n")
