import glob
import os

#### Config
project_name = config["project_name"]
terms = []
all_unique_genes = set()
for input_term in glob.glob(f"input/{project_name}/*csv"):
    term = os.path.basename(input_term).replace(".csv", "")
    terms.append(term)
    with open(input_term, "r") as f:
        genes = [l.strip() for l in f]
        config[term] = genes[1:]
        all_unique_genes.update(genes[1:])
        print(f"Added {term} with {len(genes[1:])} genes")

batches_paths = config["n_path_batches"]
batch_path_dict = dict()
genes_in_batch = int(len(all_unique_genes)/batches_paths) + 1
i = 0
b = 0

batch_path_dict[0] = []
for gene in all_unique_genes:
    i +=1
    batch_path_dict[b].append(gene)
    if i >= genes_in_batch:
        b += 1
        i = 0
        batch_path_dict[b] = []

config["path_batches"] = batch_path_dict
config["terms"] = terms
config["all_unique_genes"] = list(all_unique_genes)


#### Expected output config:
output = []
if "permutation_folder" in config:
    print("Generating permutation set, please run:")
    print(f"snakemake -s src/GraphSetup/SnakeFileDistance.smk -c 6  --configfile input/{config['permutation_folder']}.yaml")
    print("Upon completion.")
    permutations = expand("input/{permutation_folder}/{term}_set_{n}.csv",
        permutation_folder=config['permutation_folder'],
        n=range(config["permutation_N"]),
        term=config["terms"]),
    permutation_config = "input/{permutation_folder}.yaml".format(permutation_folder=config["permutation_folder"]),
    output += permutations
    output.append(permutation_config)
    include: "GenerateRandomSets/GeneratePermutSets.smk"

output.append("work/{project}/term_to_term_probability_matrix.csv".format(project = config["project_name"]))


include: "CalculatePossiblePaths/GetShortestPathNeighborhood.smk"
include: "TermNeighborhoodProbablility/TermGeneNeighborhoodProbability.smk"
include: "TermNeighborhoodProbablility/TermWithinTermNeighborhoodProbability.smk"



## Rule
rule all:
    input:
        output
