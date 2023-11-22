import glob
import os

## Config
project_name = config["project_name"]
terms = []
unique_genes = set()
for input_term in glob.glob(f"input/{project_name}/*csv"):
    term = os.path.basename(input_term).replace(".csv", "")
    terms.append(term)
    with open(input_term, "r") as f:
        genes = [l.strip() for l in f]
        config[term] = genes[1:]
        unique_genes.update(genes[1:])
        print(f"Added {term} with {len(genes[1:])} genes")

config["terms"] = terms
config["all_unique_genes"] = list(unique_genes)

include: "CalculatePossiblePaths/GetShortestPathNeighborhood.smk"
include: "TermNeighborhoodProbablility/TermGeneNeighborhoodProbability.smk"
include: "TermNeighborhoodProbablility/TermWithinTermNeighborhoodProbability.smk"


## Rule
rule all:
    input:
        "work/{project}/term_to_term_probability_matrix.csv".format(project = config["project_name"])
