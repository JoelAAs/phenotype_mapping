import multiprocessing as mp
import igraph as ig
import pandas as pd
import sys


def get_all_at_depth(input_genes, depth):
    with open(input_genes, "r") as f:
        lines = f.readlines()

    gene_depth = []
    for l in lines[1:]:
        gene = l.strip()
        for at_depth in range(1, depth+1):
            gene_depth.append([gene, at_depth])
    return gene_depth


def format_to_csv(gene_dict, depth, output):
    output.write("\t".join(map(str, range(depth))) + "\n")
    if depth == 1:
        output.write(gene_dict + "\n")
    else:
        gene = list(gene_dict.keys())[0]  # the root is always a single gene
        _format_to_line(gene, gene_dict, output, list())


def _format_to_line(gene, gene_dict, output, previous):
    previous += [gene]
    if type(gene_dict) is not dict:
        output.write("\t".join(previous) + "\n")
    else:
        for next_gene in gene_dict[gene]:
            next_previous = previous.copy()
            if type(next_gene) is dict:
                gene_name = list(next_gene.keys())[0]
                _format_to_line(gene_name, next_gene, output, next_previous)
            else:
                gene_name = next_gene
                _format_to_line(gene_name, next_gene, output, next_previous)


def get_all_neighbours_at_depth(output_folder, gene_name, depth, graph):
    print(f"getting enpoints for {gene_name} at depth {depth}")

    neighbours = _get_neighbours_at_depth(
        gene_name,
        depth,
        [],
        graph
    )
    string_id = gene_name
    output_file = open(f"{output_folder}/{string_id}_at_{depth}.csv", "w")
    format_to_csv(neighbours, depth, output_file)
    output_file.close()


def _get_neighbours_at_depth(gene_name, n, visited, graph):
    gene_vertex = graph.vs.find(gene_name)
    if n == 1:
        return gene_vertex["name"]
    else:
        visited.append(gene_vertex.index)
        current_neighbours = [
            n for n in graph.neighborhood(gene_vertex) if n not in visited
        ]
        return {
            gene_vertex["name"]: [_get_neighbours_at_depth(ni, n - 1, visited, graph) for ni in current_neighbours]}


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


def parralleled_main(input_genes, ppi_file, output_folder, n_cores=6, max_depth=6):
    gene_and_depth = get_all_at_depth(input_genes, max_depth)
    graph = setup_network(ppi_file)
    pool = mp.Pool(n_cores)
    processes = [pool.apply_async(
        get_all_neighbours_at_depth,
        args=(output_folder,
              gene_name,
              int(depth),
              graph))
        for gene_name, depth in gene_and_depth]

    for process in processes:
        process.get()
    pool.close()


if __name__ == '__main__':
    args = sys.argv[1:]
    print(args)
    input_genes = args[0]
    ppi_file = args[1]
    output_folder = args[2]
    max_depth = int(args[3])
    parralleled_main(input_genes, ppi_file, output_folder, max_depth=max_depth)
