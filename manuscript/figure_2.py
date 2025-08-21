import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
sns.set(style="whitegrid")
sc.settings.set_figure_params(dpi=300, figsize=(12, 6))  # high-res + wider figure

#load the Cortex dataset
filename = "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad"
adata = sc.read_h5ad(filename)

#filter unwanted cells
adata = adata[~adata.obs["cell_type"].isin(["NPCs", "OPCs","Microglia","Astrocytes"])].copy()

#combined group column
adata.obs["rep_cell"] = adata.obs["replicates"].astype(str) + "_" + adata.obs["cell_type"].astype(str)

#order
replicates = sorted(adata.obs["replicates"].unique(), reverse=True)
celltypes = ["Inhibitory", "Excitatory","Astrocytes"]  # your custom order
group_order = [f"{r}_{ct}" for r in replicates for ct in celltypes]

#filter only present combos
present = adata.obs["rep_cell"].unique().tolist()
group_order = [g for g in group_order if g in present]

#Plot
sc.pl.dotplot(
    adata,
    var_names=["AAV1", "AAV2", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "Retro", "Rh10", "Anc80", "SHH10"],
    groupby="rep_cell",
    use_raw=False,
    categories_order=group_order,
    dendrogram=False,
    show=True
)

sc.pl.dotplot(
    adata,
    var_names=["AAV1", "AAV2", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "Retro", "Rh10", "Anc80", "SHH10"],
    groupby="cell_type",
    use_raw=False,
    dendrogram=False,
    show=True
)


adata = sc.read_h5ad(filename)


custom_colors = {
    'Excitatory': '#0000ff',   # Blue
    'Inhibitory': '#ff6600',   # Orange
    'Astrocytes': '#008000',   # Green
    'NPCs': '#993399',         # Purple
    'OPCs': '#ff0000',         # Red
    'Microglia': '#000000',    # Black
}

# map them to your AnnData object
adata.uns['cell_type_colors'] = [
    custom_colors[c] for c in adata.obs['cell_type'].cat.categories
]

# now plot
sc.settings.set_figure_params(dpi=300, figsize=(6, 6))  # high-res + wider figure
sc.pl.umap(
    adata,
    color=["cell_type"],
    frameon=False,
    wspace=0.4,
    hspace=0.4,
    ncols=3,
    show=True
)

#use all new colors for replicates 
#red, purple, brown for replicates 1 - 3
replicate_colors = {
	'replicate_1': '#7A9F35', 
	'replicate_2': '#26596A',  
	'replicate_3': '#452F74', 
}

# map them to your AnnData object
adata.uns['replicates_colors'] = [
	replicate_colors[r] for r in adata.obs['replicates'].cat.categories
]

sc.settings.set_figure_params(dpi=300, figsize=(6, 6))  # high-res + wider figure
sc.pl.umap(
    adata,
    color=["replicates"],
    frameon=False,
    wspace=0.4,
    hspace=0.4,
    ncols=3,
    show=True
)



import numpy as np
adata = sc.read_h5ad(filename)

# Create a new log1p-transformed version of Retro
adata.obs["log1p_Retro"] = np.log1p(adata[:,"Retro"].to_df())
adata.obs["log1p_C1ql3"] = np.log1p(adata[:,"C1ql3"].to_df())

#show only excitatory and inhibitory
adata = adata[adata.obs["cell_type"].isin(["Excitatory", "Inhibitory"])].copy()

# Plot using the log-transformed values
sc.pl.umap(
    adata,
    color=["log1p_Retro",],
    frameon=False,
    wspace=0.4,
    hspace=0.4,
    ncols=3,
)

sc.pl.umap(
	adata,
	color=["log1p_C1ql3",],
	frameon=False,
	wspace=0.4,
	hspace=0.4,
	ncols=3,
)