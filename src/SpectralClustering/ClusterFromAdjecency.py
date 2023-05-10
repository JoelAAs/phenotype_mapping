from sklearn.cluster import SpectralClustering
import pandas as pd
from scipy import sparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import linalg

adj_df = pd.read_csv("data/pharmgkb/j_distance.csv", sep=" ")
adj_mat = adj_df.to_numpy()
np.fill_diagonal(adj_mat, 1)


spec_cl = SpectralClustering(
    25,
    affinity='precomputed',
    n_init=100
)
spec_cl.fit(adj_mat)
print(spec_cl.labels_)

drug_bags = dict()
with open("data/pharmgkb/bags.txt", "r") as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip().split("\t")
        i = line[0]
        drugs = line[1:]
        drug_bags[i] = drugs

def get_clusters(labels, bags):
    cluster_bags = dict()
    for i, cluster_n in enumerate(labels):
        if cluster_n in cluster_bags:
            cluster_bags[cluster_n].update(bags[str(i+1)])
        else:
            cluster_bags[cluster_n] = set(bags[str(i+1)])

    for cluster_n in cluster_bags.keys():
        drugs = cluster_bags[cluster_n]
        with open(f"data/pharmgkb/cluster/{cluster_n}.txt", "w") as w:
            for drug in drugs:
                w.write(f"{drug}\n")

    return cluster_bags


get_clusters(spec_cl.labels_, drug_bags)










### learring
sns.set_style('darkgrid', {'axes.facecolor': '.9'})
sns.set_palette(palette='deep')
sns_c = sns.color_palette(palette='deep')


adj_mat = pd.read_csv("data/pharmgkb/j_distance.csv", sep=" ")
adj_mat = adj_mat.to_numpy()
np.fill_diagonal(adj_mat, 1)


graph_laplacian = sparse.csgraph.laplacian(csgraph=adj_mat, normed=False)

# fig, ax = plt.subplots(figsize=(10, 8))
# sns.heatmap(graph_laplacian, ax=ax, cmap='viridis_r')
# ax.set(title='Graph Laplacian')


def compute_spectrum_graph_laplacian(graph_laplacian):
    """Compute eigenvalues and eigenvectors and project
    them onto the real numbers.
    """
    eigenvals, eigenvcts = linalg.eig(graph_laplacian)
    eigenvals = np.real(eigenvals)
    eigenvcts = np.real(eigenvcts)
    return eigenvals, eigenvcts


eigenvals, eigenvcts = compute_spectrum_graph_laplacian(graph_laplacian)


eigenvcts_norms = np.apply_along_axis(
  lambda v: np.linalg.norm(v, ord=2),
  axis=0,
  arr=eigenvcts
)

print('Min Norm: ' + str(eigenvcts_norms.min()))
print('Max Norm: ' + str(eigenvcts_norms.max()))

eigenvals_sorted_indices = np.argsort(eigenvals)
eigenvals_sorted = eigenvals[eigenvals_sorted_indices]

# fig, ax = plt.subplots(figsize=(10, 6))
# sns.lineplot(x=range(1, eigenvals_sorted_indices.size + 1), y=eigenvals_sorted, ax=ax)
# ax.set(title='Sorted Eigenvalues Graph Laplacian', xlabel='index', ylabel=r'$\lambda$')

index_lim = 10

fig, ax = plt.subplots(figsize=(10, 6))
sns.scatterplot(x=range(1, eigenvals_sorted_indices[: index_lim].size + 1), y=eigenvals_sorted[: index_lim], s=80, ax=ax)
sns.lineplot(x=range(1, eigenvals_sorted_indices[: index_lim].size + 1), y=eigenvals_sorted[: index_lim], alpha=0.5, ax=ax)
ax.axvline(x=6, color=sns_c[3], label='zero eigenvalues', linestyle='--')
ax.legend()
ax.set(title=f'Sorted Eigenvalues Graph Laplacian (First {index_lim})', xlabel='index', ylabel=r'$\lambda$');

zero_eigenvals_index = np.argwhere(abs(eigenvals) < 1e-5)
proj_df = pd.DataFrame(eigenvcts[:, zero_eigenvals_index.squeeze()])
proj_df.columns = ['v_' + str(c) for c in proj_df.columns]
proj_df.head()

from sklearn.cluster import KMeans
def project_and_transpose(eigenvals, eigenvcts, num_ev):
    """Select the eigenvectors corresponding to the first
    (sorted) num_ev eigenvalues as columns in a data frame.
    """
    eigenvals_sorted_indices = np.argsort(eigenvals)
    indices = eigenvals_sorted_indices[: num_ev]

    proj_df = pd.DataFrame(eigenvcts[:, indices.squeeze()])
    proj_df.columns = ['v_' + str(c) for c in proj_df.columns]
    return proj_df

inertias = []

k_candidates = range(1, 6)

for k in k_candidates:
    k_means = KMeans(random_state=42, n_clusters=k)
    k_means.fit(proj_df)
    inertias.append(k_means.inertia_)

    fig, ax = plt.subplots(figsize=(10, 6))
sns.scatterplot(x=k_candidates, y = inertias, s=80, ax=ax)
sns.lineplot(x=k_candidates, y = inertias, alpha=0.5, ax=ax)
ax.set(title='Inertia K-Means', ylabel='inertia', xlabel='k')

k = 6

def run_k_means(df, n_clusters):
    """K-means clustering."""
    k_means = KMeans(random_state=25, n_clusters=n_clusters)
    k_means.fit(df)
    cluster = k_means.predict(df)
    return cluster

cluster = run_k_means(proj_df, n_clusters=6)
