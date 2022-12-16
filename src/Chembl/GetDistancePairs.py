import igraph as ig
import pandas as pd
import ensembl_rest
from itertools import combinations
import requests
from ensembl_rest import HTTPError
import numpy as np
import multiprocessing as mp


def _get_depth_of_pair(pair, graph):
    depth = len(graph.get_all_shortest_paths(pair[0], pair[1])[0])
    return depth


class GetGeneDistanceAtN:
    def __init__(self, ppi_network_file, output_folder, output_names, *args):
        self.graph = self.setup_network(ppi_network_file)
        self.output_folder = output_folder
        self.output_names = output_names
        self.args = args
        self.gene_names = []
        self.string_ids = []
        self.network_index = []
        self.name_df = []
        self.gene_pairs = []
        self.distance_pairs = []

    def calculated_distance_pairs(self):
        self.gene_names = self.get_all_human_genes_set(self.args)
        print(f"{len(self.gene_names)} unique genes loaded")
        self.string_ids = self.get_stringdb_identifiers()
        self.network_index = self.get_network_index()
        self.name_df = self.write_name_df(self.output_names)
        self.gene_pairs = [(a, b) for a, b in combinations(self.network_index, 2)]
        self.distance_pairs = self.set_depth_of_pairs()

    def load_previous_distances(self, distance_file, name_df):
        distance_pairs = []
        with open(distance_file, "r") as f:
            for line in f:
                print(line)
                pair_1, pair_2, distance = line.strip().split("\t")
                distance_pairs.append(((int(pair_1), int(pair_2)), int(distance)))
        self.distance_pairs = distance_pairs
        self.name_df = pd.read_csv(name_df, sep="\t")

    """
    Gets all indexes for nodes based on ensembl ids
    """

    def get_network_index(self):
        network_indexes = []
        for string_id in self.string_ids:
            index = self.graph.vs.find(string_id).index
            network_indexes.append(index)
        return network_indexes

    def get_stringdb_identifiers(self):
        string_ids = []
        for gene in self.gene_names:
            print(f"looking up string id for {gene}")
            url = "https://string-db.org/api/json/get_string_ids?identifiers={gene}&species=9606"
            response = requests.get(url.format(gene=gene))
            if response.ok:
                data = response.json()
                string_ids.append(data[0]["stringId"])
            else:
                raise HTTPError(response.text + f" for {gene}")

        return string_ids

    """
    Queries Ensembl for ids based on gene names from Chembl
    """

    def get_ensembl_ids(self):
        ensembl_ids = []
        for gene in self.gene_names:
            try:
                print(f"looking up ensembl id for {gene}")
                response = ensembl_rest.symbol_lookup(
                    species="homo sapiens",
                    symbol=gene
                )
                ensembl_ids.append(response["id"])  # ensembl id
            except HTTPError:
                synonym_url = "https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene}?"
                response = requests.get(
                    synonym_url.format(gene=gene),
                    headers={"Content-Type": "application/json"}
                )
                if response.ok:
                    data = response.json()
                else:
                    raise HTTPError(f"Response for {gene}: {response.status}")
                if len(data) > 1:
                    raise HTTPError(f"Multiple hits for {gene}")
                elif data[0]["type"] != "gene":
                    raise HTTPError(f"{gene} didn't get a hit among genes")
                else:
                    ensembl_ids.append(data[0]["id"])  # ensembl id

        return ensembl_ids

    """
    Gets all gene names from targets
    """

    @staticmethod
    def get_all_human_genes_set(df_tuple):
        gene_set = set()
        for arg in df_tuple:
            gene_interaction_df = pd.read_csv(arg, sep="\t")
            gene_set.update(
                {gene for gene in gene_interaction_df["gene"] if gene != "Non-human"}
            )

        return list(gene_set)

    # We assume node/node/weight columns
    @staticmethod
    def setup_network(ppi_network_file):
        edge_list_df = pd.read_csv(ppi_network_file, sep=" ")
        edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
        edge_list_df = edge_list_df[edge_list_df["combined_score"] > 0.6]
        tuple_list = [tuple(x) for x in edge_list_df.values]
        graph = ig.Graph.TupleList(
            tuple_list,
            directed=False,
            weights=True)
        return graph

    def set_depth_of_pairs(self):
        with mp.Pool(mp.cpu_count() - 1) as pool:
            star_args = [(pair, self.graph) for pair in self.gene_pairs]
            distances = pool.starmap(
                _get_depth_of_pair, star_args
            )

        distance_pairs = list(zip(self.gene_pairs, distances))

        return distance_pairs

    def get_all_at_depth(self):
        gene_per_depth = set()
        for pair, depth in self.distance_pairs:
            gene_a, gene_b = pair
            gene_per_depth.update({(gene_a, depth), (gene_b, depth)})
        return list(gene_per_depth)

    def get_all_neighbours_at_depth(self):
        gene_and_depth = self.get_all_at_depth()
        for gene_index, depth in gene_and_depth:
            print(f"getting enpoints for {gene_index} at depth {depth}")
            neighbours = self._get_neighbours_at_depth(
                gene_index,
                depth,
                []
            )
            string_id = self.graph.vs.find(gene_index)["name"]
            output_file = open(f"{self.output_folder}/{string_id}_at_{depth}.csv", "w")
            self.format_to_csv(neighbours, depth, output_file)
            output_file.close()

    def _get_neighbours_at_depth(self, gene_index, n, visited):
        gene_vertex = self.graph.vs.find(gene_index)
        if n == 1:
            return gene_vertex["name"]
        else:
            visited.append(gene_vertex.index)
            current_neighbours = [
                n for n in self.graph.neighborhood(gene_vertex) if n not in visited
            ]
            return {
                gene_vertex["name"]: [self._get_neighbours_at_depth(ni, n - 1, visited) for ni in current_neighbours]}

    def write_name_df(self, output_names, write=True):
        df = pd.DataFrame(
            {"gene_name": self.gene_names,
             "string_id": self.string_ids,
             "network_index": self.network_index
             }
        )
        if write:
            df.to_csv(output_names, index=False, sep="\t")
        return df

    def format_to_csv(self, gene_dict, depth, output):
        output.write("\t".join(map(str, range(depth))) + "\n")
        gene = list(gene_dict.keys())[0]  # the root is always a single gene
        self._format_to_line(gene, gene_dict, output, list())

    def _format_to_line(self, gene, gene_dict, output, previous):
        previous += [gene]
        if type(gene_dict) is not dict:
            output.write("\t".join(previous) + "\n")
        else:
            for next_gene in gene_dict[gene]:
                next_previous = previous.copy()
                if type(next_gene) is dict:
                    gene_name = list(next_gene.keys())[0]
                    self._format_to_line(gene_name, next_gene, output, next_previous)
                else:
                    gene_name = next_gene
                    self._format_to_line(gene_name, next_gene, output, next_previous)

    def write_gene_distances(self, output_distance_file):
        with open(output_distance_file, "w") as w:
            for pair, distance in self.distance_pairs:
                w.write(f"{pair[0]}\t{pair[1]}\t{distance}\n")



if __name__ == '__main__':
    import glob

    args = glob.glob("work/predicted_interaction/chembl/annotated/*csv")
    ppi_network_file = "data/9606.protein.links.v11.5.txt"
    output_folder = "test"
    output_names = "work/chembl/pairs/name_conversion.csv"
    previous_distances = "test/pair_distance.csv"

    Gt = GetGeneDistanceAtN(ppi_network_file, output_folder, output_names, *args)
    Gt.calculated_distance_pairs()
    Gt.write_gene_distances("work/chembl/pairs/distance_pairs.csv")
