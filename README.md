## Introduction
This analysis is used for abstraction of distances of sets of genes based on their positions in a interaction network.
The method were developed for clustering sets of genes associated with adverse drug reactions and drugs.

Given the parameters set (set in `input/myproject.yaml`) it will cluster the different sets of terms (read drugs/ADRs) based on the shared interaction paths in the interaction network.
For each cluster, the proteins are ranked on importance for the clusters and a subsequent KEGG enrichment analysis is performed.

### Dependencies and versions
Code run on the following packages:
- `Snakemake 7.18.2`
- `Python 3.9.7`
  - `bzip2 1.0.8`
  - `igraph 0.10.2`
  - `pandas 1.4.0`
  - `scipy 1.9.3`
  - `scikit-learn 1.0.2`
  - `networkx 3.0`
- `R 4.3.2`
  - `ggplot2 3.4.4`
  - `ggrepel 0.9.4`
  - `dplyr 1.1.4`
  - `gridExtra 2.3`
  - `clusterProfiler 4.10.1`

### Data and input
#### Data
Please download human protein-protein interactions from StringDB.
Place the links file in `data`.

#### Input:
In the folder `input` write a yaml configfile for your project and set parameters as:
```
project_name: "example_input"                 # write your name of the project, this is also the name of the folder containing your term-gene csvs
ppi_file: "data/9606.protein.links.v11.5.txt" # path to the stringDB interaction file
combined_score_threshold: 0.7                 # Sets the threshold for the combined_score in the StringDB
max_depth: 4                                  # The maximum length of interaction paths considerd 
n_path_batches: 1                             # set number of batches for interaction paths estimation (default 6 cores per batch) 
n_permutation_combinations: 10                # Number of equivelent sets per term, from with baseline probabilities ar calculated

n_clusters: 3                             # Number clusters estimated (will crash if there are less terms than clusters) 
permutation_folder: "example_permutation" # The project for the calulations of baseline probabilities
permutation_N: 100                         # Number of combinations of n_permutation_combinations will be estimated 
```

Under `input` create a folder with the same name as `project_name` in the configfile.
Place the list in that file with the genes expressed as StringDB identifiers, example `input/my_project/my_drug1.csv`:
```
gene
9606.ENSP00000295897
9606.ENSP00000259396
...
```

## Exceution
### GraphSetup
First run:
```
snakemake -s src/GraphSetup/SnakeFileDistance.smk -c 6  --configfile input/my_config.yaml
```

This will create `work/my_project/term_to_term_probability_matrix.csv` containing the proximities between all input-terms.
In addition, a folder containing all equivalent random sets for terms will are located in `input/my_permutations/` and a yaml file `input/my_permutations.yaml`.

Run:
```
snakemake -s src/GraphSetup/SnakeFileDistance.smk -c 6  --configfile input/my_permutations.yaml
```
This will calculate baseline probabilities.

### GraphAnalysis
Running:
```
snakemake -s src/GraphAnalysis/SnakeFileanalysis.smk -c 6  --configfile input/my_config.yaml
```
Will cluster the terms and save the clusters under `work/my_config/clustering/SCnorm_{n_clusters}.csv"`.
Subsequent files can be found:
```
work/{project}/group-quant_{n_clusters}/{cluster}_quant.csv  # Ranking of genes in cluster {cluster}
work/{project}/group-quant_{n_clusters}/{cluster}_top.csv # Top 500 ranking genes (can be changes in rule get_top_values)
work/{project}/group-quant_{n_clusters}/{cluster}_connected_components_enrichment/Component_{component}.csv  # KEGG enrichment results of a component in {cluster}
work/{project}/group-quant_{n_clusters}/{cluster}_connected_components_enrichment/Component_{component}_KEGG.png"  # Graph of enrichment results
```

### Replication
For replication of results please follow instructions of listed in https://github.com/JoelAAs/phenotype_mapping/tree/replication
