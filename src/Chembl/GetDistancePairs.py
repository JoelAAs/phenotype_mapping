import igraph as ig
import pandas as pd
import ensembl_rest
import requests
from ensembl_rest import HTTPError
import multiprocessing as mp
import os


def _write_part_depth(pairs, graph, pro_n, name_dict):
    i = 0
    with open(f"tmp/part_{pro_n}", "w") as w:
        for pair in pairs:
            depth = len(graph.get_all_shortest_paths(pair[0], pair[1])[0])
            i += 1
            if i % 100 == 0:
                print(f"Finished {round(i / len(pairs) * 100)} % of process {pro_n}")
            w.write(f"{name_dict[pair[0]]}\t{name_dict[pair[1]]}\t{depth}\n")


class GetGeneDistanceAtN:
    def __init__(self, ppi_network_file, output_folder, output_names, combination_file, missing_filename):
        self.graph = self.setup_network(ppi_network_file)
        self.output_folder = output_folder
        self.output_names = output_names
        self.combination_file = combination_file
        self.named_gene_pairs = None
        self.gene_names = []
        self.string_ids = []
        self.translation_dict_gene = None
        self.translation_dict_index = None
        self.translation_dict_string = None
        self.network_index = []
        self.name_df = []
        self.gene_pairs = []
        self.distance_pairs = []
        self.missing = missing_filename

    def calculated_distance_pairs(self):
        self.named_gene_pairs = self.get_all_human_genes_set(self.combination_file)
        self.gene_names = list(set([gene_name for pair in self.named_gene_pairs for gene_name in pair]))
        print(f"{len(self.gene_names)} unique genes loaded")
        self.string_ids, self.translation_dict_gene = self.get_stringdb_identifiers()
        self.network_index, self.translation_dict_index, self.translation_dict_string = self.get_network_index()
        self.remove_failed_names()
        self.name_df = self.write_name_df(self.output_names)
        t = lambda x: self.translation_dict_index[self.translation_dict_gene[x]]
        self.gene_pairs = [(t(a), t(b)) for a, b in self.named_gene_pairs]
        self.distance_pairs = self.set_depth_of_pairs()

    def remove_failed_names(self):
        failed = pd.read_csv(self.missing, sep="\t")
        pairs = []
        for a, b in self.named_gene_pairs:
            if not (a in failed.name or b in failed.name):
                try:
                    _ = self.translation_dict_index[self.translation_dict_gene[a]]
                    _ = self.translation_dict_index[self.translation_dict_gene[b]]
                    pairs.append((a, b))
                except KeyError:
                    pass
        self.named_gene_pairs = pairs

    def load_previous_distances(self, distance_file, name_df):
        distance_pairs = []
        self.name_df = pd.read_csv(name_df, sep="\t")
        translation_dict = dict()
        for i, gene_id, string_id, index_id in self.name_df.itertuples():
            translation_dict.update({string_id: int(index_id)})
        self.translation_dict_index = translation_dict
        with open(distance_file, "r") as f:
            for line in f:
                pair_1, pair_2, distance = line.strip().split("\t")
                distance_pairs.append((
                    (self.translation_dict_index[pair_1],
                     self.translation_dict_index[pair_2]), int(distance)))
        print(distance_pairs)
        self.distance_pairs = distance_pairs

    """
    Gets all indexes for nodes based on ensembl ids
    """

    def get_network_index(self):
        network_indexes = []
        translation_dict = dict()
        translation_string = dict()
        missing = open(self.missing, "a")
        for string_id in self.string_ids:
            try:
                index = self.graph.vs.find(string_id).index
                network_indexes.append(index)
                translation_dict.update({string_id: index})
                translation_string.update({index: string_id})
            except ValueError:
                missing.write(f"{string_id}\tvertex missing\n")

        missing.close()
        return network_indexes, translation_dict, translation_string

    def get_stringdb_identifiers(self):
        string_ids = []
        translation_dict = dict()
        url = "https://string-db.org/api/json/get_string_ids?identifiers={gene}&species=9606"
        missing = open(self.missing, "w")
        missing.write("name\treason\n")
        print(f"there are {len(self.gene_names)} unique genes")
        for gene in self.gene_names:
            response = requests.get(url.format(gene=gene))
            print(f"looking up string id for {gene}")
            if response.ok:
                data = response.json()
                if data:
                    sid = data[0]["stringId"]
                    string_ids.append(sid)
                    translation_dict.update({gene: sid})
                else:
                    missing.write(f"{gene}\tnot in stringdb\n")
            else:
                raise HTTPError(response.text + f" for {gene}")
        missing.close()

        return string_ids, translation_dict

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
    def get_all_human_genes_set(combination_file):
        df_combinations = pd.read_csv(combination_file, sep="\t")
        gene_interactions = [(a, b) for i, a, b in df_combinations.itertuples()]
        return gene_interactions

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
        os.mkdir("tmp")
        n_cores = mp.cpu_count()

        def _binit(inlist, steps, step_len, current_step, outlist):
            if current_step == steps:
                return outlist + [inlist]
            outlist = outlist + [inlist[0:step_len]]
            return _binit(inlist[step_len:], steps, step_len, current_step + 1, outlist)

        parts = 100
        binned_pairs = _binit(self.gene_pairs, parts, int(len(self.gene_pairs) / parts), 0, [])

        pool = mp.Pool(n_cores)
        processes = [pool.apply_async(
            _write_part_depth,
            args=(part,
                  self.graph,
                  i,
                  self.translation_dict_string))
            for i, part in enumerate(binned_pairs)]

        for process in processes:
            process.get()
        pool.close()

        distance_pairs = []
        for i in range(parts + 1):
            with open(f"tmp/part_{i}", "r") as f:
                for line in f:
                    gene_from, gene_to, depth = line.strip().split("\t")
                    distance_pairs.append(((gene_from, gene_to), int(depth)))

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
        print("writing name file")
        rows = []
        for gene in self.gene_names:
            try:
                rows.append({
                    "gene_name": gene,
                    "string_id": self.translation_dict_gene[gene],
                    "network_index": self.translation_dict_index[self.translation_dict_gene[gene]]
                })

            except KeyError:
                pass

        df = pd.DataFrame(rows)
        if write:
            df.to_csv(output_names, index=False, sep="\t")
        return df

    def format_to_csv(self, gene_dict, depth, output):
        output.write("\t".join(map(str, range(depth))) + "\n")
        if depth == 1:
            output.write(gene_dict + "\n")
        else:
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
        print("writing distance file")
        with open(output_distance_file, "w") as w:
            for pair, distance in self.distance_pairs:
                w.write(f"{pair[0]}\t{pair[1]}\t{distance}\n")
        print("done writing distance file")
