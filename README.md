## Introduction
In order to replicate results published please execute the following steps.

### Dependencies and versions
Code ran on the following packages:
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

### Data
Please download human protein-protein interactions from StringDB. This analysis was done on version 11.5.
Place the links file in `data`. If using another version change the path to the file among the snakemake rules or rename it to `data/9606.protein.links.v11.5.txt`
## Replication

### Calculate interaction based edge affinity for ADR and drug networks
The analysis for the interaction based edge weights are defined by code under `src/GraphSetup`.

![Rulegraph of GraphSetup](GraphSetup.png)

First, calculate the drug-drug probabilities:
```
snakemake -s src/GraphSetup/SnakeFileDistance.smk -c 1 --configfile input/full-drugbank.yaml
```

Secondly, calculate the ADR-ADR probabilities:
```
snakemake -s src/GraphSetup/SnakeFileDistance.smk -c 6 --configfile input/HPO-pruned.yaml
```
This step will take ~4-5 min to calculate the DAG. The whole analysis workflow should take around ~24 h on six core.

Join data into the bipartite drug-ADR network by combining the edges as:
```
mkdir work/MedAdr
cat work/full-drugbank/term_to_term_probability_matrix.csv > work/MedAdr/term_to_term_probability_matrix.csv
tail -n +2 work/HPO-pruned/term_to_term_probability_matrix.csv >> work/MedAdr/term_to_term_probability_matrix.csv
tail -n +2 data/hpo_med_directed.csv >> work/MedAdr/term_to_term_probability_matrix.csv
```

Generate the probability per the gene given random sets per drug for CDF calculations later:
```
snakemake -s src/GraphSetup/SnakeFileDistance.smk -c 6 --configfile input/full-drugbank-permut.yaml
```
The random sets for each drug were generated by `src/GraphAnalysis/CandidateGeneAnalysis/GeneratePermutSets.smk`.
This step will take a long time and if computational resources are available I recommend using them.
Please set the `n_path_batches` parameter to the desired number of batches to run in parallel.
Observe that `src/GraphSetup/CalculateShortestPaths/possible_paths_from_a.py` is set to a default of 6 cores, that is 6 cores per batch.

### Phenotype clustering and candidate genes
#### Phenotype clustering 
First we calculate the `CDF(Zd)` affinity, the `Sab` affinity and normalised probabilities as:
```
snakemake -s src/GraphAnalysis/SnakeFileanalysis.smk -c 6  \
 work/full-drugbank-benchmark/term_to_term_probability_matrix_zscore.csv \
 work/HPO-pruned-benchmark/term_to_term_probability_matrix_zscore.csv \
 work/full-drugbank/term_to_term_probability_matrix_norm.csv
```

First we obtain the clusters from `full-drugbank` and `MedAdr`:
```
snakemake -s src/GraphAnalysis/SnakeFileanalysis.smk -c 1  work/full-drugbank/clustering/SCnorm_3.csv  work/MedAdr/clustering/SCnorm_3.csv
```

Similar to earlier we join the drug and ADR affinities:
```
mkdir work/MedAdr-benchmark
cat work/full-drugbank-benchmark/term_to_term_probability_matrix_zscore.csv > work/MedAdr-benchmark/term_to_term_probability_matrix_zscore.csv
cat work/full-drugbank-benchmark/term_to_term_probability_matrix_sab.csv > work/MedAdr-benchmark/term_to_term_probability_matrix_sab.csv
cat work/full-drugbank/term_to_term_probability_matrix_norm.csv > work/MedAdr/term_to_term_probability_matrix_norm.csv
tail -n +2 work/HPO-pruned-benchmark/term_to_term_probability_matrix_zscore.csv >> work/MedAdr-benchmark/term_to_term_probability_matrix_zscore.csv
tail -n +2 work/HPO-pruned-benchmark/term_to_term_probability_matrix_sab.csv >> work/MedAdr-benchmark/term_to_term_probability_matrix_sab.csv
tail -n +2 work/HPO-pruned/term_to_term_probability_matrix_norm.csv >> work/MedAdr/term_to_term_probability_matrix_norm.csv
tail -n +2 data/hpo_med_directed.csv >> work/MedAdr-benchmark/term_to_term_probability_matrix_zscore.csv
tail -n +2 data/hpo_med_directed.csv >> work/MedAdr-benchmark/term_to_term_probability_matrix_sab.csv
tail -n +2 data/hpo_med_directed.csv >> work/MedAdr/term_to_term_probability_matrix_norm.csv
```

Running:
```
snakemake -s src/GraphAnalysis/SnakeFileanalysis.smk -c 6 work/meta_plots/Sensitivity.png
```
Gives us the sensitivity of each method and the corresponding graph at `work/meta_plots/Sensitivity.png`.
The specific sensitivity files are found at:
```
work/full-drugbank/clustering/metrics/sensitivity_inbew.csv
work/full-drugbank/clustering/metrics/sensitivity_norm.csv
work/full-drugbank-benchmark/clustering/metrics/sensitivity_sab.csv
work/full-drugbank-benchmark/clustering/metrics/sensitivity_zscore.csv
work/MedAdr/clustering/metrics/sensitivity_inbew.csv
work/MedAdr/clustering/metrics/sensitivity_norm.csv
work/MedAdr-benchmark/clustering/metrics/sensitivity_sab.csv
work/MedAdr-benchmark/clustering/metrics/sensitivity_zscore.csv
work/MedAdr/clustering/metrics/hpo_sensitivity_inbew.csv
work/MedAdr/clustering/metrics/hpo_sensitivity_norm.csv
work/MedAdr-benchmark/clustering/metrics/hpo_sensitivity_sab.csv
work/MedAdr-benchmark/clustering/metrics/hpo_sensitivity_zscore.csv
```

#### Candidate genes
To generate the `CDF(InBEW)` we first need to calculate the probabilities of drawing each gene from the clusters.
There are previous probabilities in the repository but to replicate these run:
```
snakemake -s src/GraphAnalysis/SnakeFileanalysis.smk -c 6 work/full-drugbank/candidate_genes/probabilities_3
```

Match the three files (`cluser_0.csv` to `cluster_2.csv` in `work/full-drugbank/candidate_genes/probabilities_3`) against the drug groups in each cluster (listed in `work/full-drugbank/clustering/SCnorm_3.csv`).
The drug groups are specified in `data/Node_groups.csv`.

Copy the correct cluster to each named files in `data/full_drugbank_gene_probabilites` as:
```
cp work/full-drugbank/candidate_genes/probabilities_3/cluster_{N}.csv  data/full_drugbank_gene_probabilites/{group}_gene_probabilites.csv
```

For all three clusters.

Now run (this step is quite memory intensive, if OOM decrease `n_permutation_combinations` specified in `SnakeFileanalysis.smk`):
```
snakemake -s src/GraphAnalysis/SnakeFileanalysis.smk -c 1
``` 

Enrichment analysis for the large connected components och each drug group can be found under `work/full-drugbank-benchmark/group-quant/{group}_connected_components_enrichment`.

For DIAMOnD graphs please run `paper/manual_figures/diamond_vs_CDF.R`.
