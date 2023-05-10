import requests
import pandas as pd
import igraph as ig


def get_string_ids(stringdb_name):
    translate_dict = dict()
    with open(stringdb_name, "r") as f:
        head = True
        for line in f:
            if head:
                head = False
            else:
                line = line.strip().split("\t")
                translate_dict.update({line[1]: line[0]})
    return translate_dict


def check_if_vertex_passes_score(graph, vertex_name):
    return len(graph.vs.select(name=vertex_name)) != 0


def setup_network(ppi_network_file, min_score=0.6):
    edge_list_df = pd.read_csv(ppi_network_file, sep=" ")
    edge_list_df["combined_score"] = edge_list_df["combined_score"] / 1000
    edge_list_df = edge_list_df[edge_list_df["combined_score"] > min_score]
    tuple_list = [tuple(x) for x in edge_list_df.values]
    ppi_graph = ig.Graph.TupleList(
        tuple_list,
        directed=False,
        weights=True)
    return ppi_graph


def query_hpo(hpofile, output_location, stringdb_name, stringdb_network):
    error_log = open(f"{output_location}/error.log", "w")
    hpo_df = pd.read_csv(hpofile, sep="\t")
    hpoids = hpo_df.HPO.unique()
    hpoids = hpoids[hpoids == hpoids]
    translation_dict = get_string_ids(stringdb_name)
    graph = setup_network(stringdb_network)

    for hpoid in hpoids:
        hpoid_file_name = hpoid.replace(":", "-")
        with open(f"{output_location}/{hpoid_file_name}.csv", "w") as w:
            w.write("gene\n")
            print(f"Querying {hpoid}")
            if "ORPHA" in hpoid:
                api = "https://hpo.jax.org/api/hpo/disease/{hpoid}?max=-1"
                key = "geneAssoc"
            else:
                api = "https://hpo.jax.org/api/hpo/term/{hpoid}/genes?max=-1"
                key = "genes"
            query = api.format(hpoid=hpoid)
            response = requests.get(query, verify="data/jax-org-chain.pem")
            if response.ok:
                data = response.json()
                for gene in data[key]:
                    gene_name = gene["geneSymbol"]
                    if gene_name in translation_dict:
                        string_id = translation_dict[gene_name]
                        if check_if_vertex_passes_score(graph, string_id):
                            w.write(f"{string_id}\n")
                        else:
                            error_log.write(f"{gene_name} not high enough interaction score\n")
                    else:
                        error_log.write(f"{gene_name} not in stringdb\n")
            else:
                error_log.write(f"{hpoid} Server: {response.status_code}\n")
                error_log.close()


if __name__ == '__main__':
    hpofile = "data/micromedex/all_unique_terms.csv"
    output_location = "input/HPO"
    stringdb_name = "data/stringdb/9606.protein.info.v11.5.txt"
    stringdb_network = "data/9606.protein.links.v11.5.txt"
    query_hpo(hpofile, output_location, stringdb_name, stringdb_network)