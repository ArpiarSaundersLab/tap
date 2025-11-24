import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tap as t
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
sns.set(style="whitegrid")
sc.settings.set_figure_params(dpi=300)

#opens the Cortex dataset
filename = "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad"
adata = sc.read_h5ad(filename)

#selects only replicate 1
adata = adata[adata.obs["replicates"]=="replicate_1"].copy()

#generates a TAP of the Cortex dataset \\ replicate 1 \\ excitatory and inhibitory only
parameters = {
	"adataObject": adata,
	"name": "Cortex Single",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"exclude" : ["Microglia","NPCs","OPCs","Astrocytes","Oligodendrocytes"],
	"outputPath": "runs/",
	"outputName": "Cortex_Neurons_Only_Single.html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler",
	"replicate_mode": False,
	"excludeMarkers" : True,
}
results = t.TAP(**parameters)


#here, we get a fresh data copy
filename = "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad"
adata = sc.read_h5ad(filename)

#filters unwanted cells
adata = adata[~adata.obs["cell_type"].isin(["NPCs", "OPCs","Microglia","Astrocytes","Oligodendrocytes"])].copy()

#remove replicate 4
adata = adata[adata.obs["replicates"]!="replicate_4"].copy()

adata.obs.replicates.value_counts()

#generates a TAP of the Cortex dataset all replicates \\ excitatory and inhibitory only
parameters = {
	"adataObject": adata,
	"name": "Cortex Reps 1-3",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	"exclude" : ["Microglia","NPCs","OPCs","Astrocytes","Oligodendrocytes"],
	"outputPath": "runs/",
	"outputName": "Cortex_Neurons_Only_All.html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 5,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler",
	"replicate_mode": True,
	"excludeMarkers" : True,
}
results = t.TAP(**parameters)


#fresh copy of data
filename = "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad"
adata = sc.read_h5ad(filename)

#selects only replicate 1
adata = adata[adata.obs["replicates"]=="replicate_1"].copy()

#filters unwanted cells
adata = adata[~adata.obs["cell_type"].isin(["NPCs", "OPCs","Microglia","Astrocytes"])].copy()

#combines the group column
adata.obs["rep_cell"] = adata.obs["replicates"].astype(str) + "_" + adata.obs["cell_type"].astype(str)

#orders the groups
replicates = sorted(adata.obs["replicates"].unique(), reverse=True)
celltypes = ["Inhibitory", "Excitatory","Astrocytes"] 
group_order = [f"{r}_{ct}" for r in replicates for ct in celltypes]

#filter only present combos
present = adata.obs["rep_cell"].unique().tolist()
group_order = [g for g in group_order if g in present]

#plots with grouped cell types with replicates
sc.pl.dotplot(
    adata,
    var_names=["AAV1", "AAV2", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "Retro", "Rh10", "Anc80", "SHH10"],
    groupby="rep_cell",
    use_raw=False,
    categories_order=group_order,
    dendrogram=False,
    show=True
)

#grouped cell types for all replicates
sc.pl.dotplot(
    adata,
    var_names=["AAV1", "AAV2", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "Retro", "Rh10", "Anc80", "SHH10"],
    groupby="cell_type",
    use_raw=False,
    dendrogram=False,
    show=True
)


#reloads the data
adata = sc.read_h5ad(filename)
adata = adata[adata.obs["replicates"]!="replicate_4"].copy()

#define a custom color palette
custom_colors = {
    'Excitatory': '#0000ff',   #Blue
    'Inhibitory': '#ff6600',   #Orange
    'Astrocytes': '#008000',   #Green
    'NPCs': '#993399',         #Purple
    'OPCs': '#ff0000',         #Red
    'Microglia': '#000000',    #Black
}

#map to AnnData object
adata.uns['cell_type_colors'] = [
    custom_colors[c] for c in adata.obs['cell_type'].cat.categories
]

#cell type umap with custom colors
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

#uses all new colors for replicates 
replicate_colors = {
	'replicate_1': '#7A9F35', 
	'replicate_2': '#26596A',  
	'replicate_3': '#452F74', 
}


#maps to AnnData object
adata.uns['replicates_colors'] = [
	replicate_colors[r] for r in adata.obs['replicates'].cat.categories
]

#replicates umap
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

#c1ql3 umap
sc.settings.set_figure_params(dpi=300, 
figsize=(6, 6), 
fontsize=20) 

sc.pl.umap(
	adata,
	color=["C1ql3"],
	frameon=False,
	wspace=0.4,
	hspace=0.4,
	ncols=2,
	show=True,
	use_raw=False,
	s=10,
)

tap_colors = [(0.0, 'blue'), (0.1, 'white'), (1.0, 'red')]
tap_cmap = LinearSegmentedColormap.from_list('tap_whitebluered', tap_colors, N=256)

#retro umap
sc.pl.umap(
	adata,
	color=["Retro"],
	frameon=False,
	wspace=0.4,
	hspace=0.4,
	ncols=2,
	show=True,
	use_raw=False,
	s=10,
	cmap=tap_cmap
)


#reloads the data 
adata = sc.read_h5ad(filename)
adata = adata[adata.obs["replicates"]!="replicate_4"].copy()

#filters excitatory neurons
excitatory_data = adata[adata.obs["cell_type"]=="Excitatory"].copy()

#removes outliers using IQR method for C1ql3 expression
c1ql3_expr = excitatory_data[:, "C1ql3"].X.toarray().flatten()
Q1 = np.percentile(c1ql3_expr, 25)
Q3 = np.percentile(c1ql3_expr, 75)
IQR = Q3 - Q1
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

#keeps only cells within the outlier bounds
outlier_mask = (c1ql3_expr >= lower_bound) & (c1ql3_expr <= upper_bound)
excitatory_filtered = excitatory_data[outlier_mask].copy()

#add replicate colors to the filtered data
excitatory_filtered.uns['replicates_colors'] = [
	replicate_colors[r] for r in excitatory_filtered.obs['replicates'].cat.categories
]

#voilin plots of c1ql3 in excitatory neurons for each replicate (outliers removed)
sc.settings.set_figure_params(dpi=300, 
figsize=(6, 4), 
fontsize=20) 
sc.pl.violin(
	excitatory_filtered,
	keys="C1ql3",
	groupby="replicates",
	stripplot=False,
	show=True,
	use_raw=True
)

#removes outliers using IQR method for Retro expression
retro_expr = excitatory_data[:, "Retro"].X.toarray().flatten()
Q1_retro = np.percentile(retro_expr, 25)
Q3_retro = np.percentile(retro_expr, 75)
IQR_retro = Q3_retro - Q1_retro
lower_bound_retro = Q1_retro - 1.5 * IQR_retro
upper_bound_retro = Q3_retro + 1.5 * IQR_retro

#keeps only cells within the outlier bounds for Retro
retro_outlier_mask = (retro_expr >= lower_bound_retro) & (retro_expr <= upper_bound_retro)
excitatory_retro_filtered = excitatory_data[retro_outlier_mask].copy()

#add replicate colors to the filtered data
excitatory_retro_filtered.uns['replicates_colors'] = [
	replicate_colors[r] for r in excitatory_retro_filtered.obs['replicates'].cat.categories
]


#violin plot for retro (outliers removed)
sc.settings.set_figure_params(dpi=300, 
figsize=(6, 4), 
fontsize=20)
sc.pl.violin(
	excitatory_retro_filtered,
	keys="Retro",
	groupby="replicates",
	stripplot=False,
	show=True,
	use_raw=True
)

###########################################################################
# Dotplot of all cell types in replicates 1,2,3
###########################################################################
#reloads the data
adata = sc.read_h5ad(filename)
adata = adata[adata.obs["replicates"]!="replicate_4"].copy()

#dotplot of all cell types in replicates 1,2,3
sc.pl.dotplot(
	adata,
	var_names=["AAV1", "AAV2", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "Retro", "Rh10", "Anc80", "SHH10"],
	groupby="cell_type",
	use_raw=False,
	dendrogram=False,
	show=False,
	swap_axes=True
)
plt.yticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()




###########################################################################
# SUPPLEMENTARY FIGURES
###########################################################################

#reloads the data
adata = sc.read_h5ad(filename)
adata1 = adata[adata.obs["replicates"]!="replicate_1"].copy()
adata2 = adata[adata.obs["replicates"]!="replicate_2"].copy()
adata3 = adata[adata.obs["replicates"]!="replicate_3"].copy()

#stacked_violin for each replicate
sc.pl.stacked_violin(
	adata1,
	var_names=["AAV1", "AAV2", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "Retro", "Rh10", "Anc80", "SHH10"],
	groupby="cell_type",
	use_raw=False,
	dendrogram=False,
	swap_axes=True,
	show=True
)

sc.pl.stacked_violin(
	adata2,
	var_names=["AAV1", "AAV2", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "Retro", "Rh10", "Anc80", "SHH10"],
	groupby="cell_type",
	use_raw=False,
	dendrogram=False,
	swap_axes=True,
	show=True
)

sc.pl.stacked_violin(
	adata3,
	var_names=["AAV1", "AAV2", "AAV5", "AAV6", "AAV7", "AAV8", "AAV9", "Retro", "Rh10", "Anc80", "SHH10"],
	groupby="cell_type",
	use_raw=False,
	dendrogram=False,
	swap_axes=True,
	show=True
)

###########################################################################
# Venn  diagram of  retro excitatory features
###########################################################################

#imports
import pandas as pd

#io
feature_importance_df=pd.read_csv("runs/feature_importances_replicates.csv")

#drop stray index cols
df=feature_importance_df.loc[:, ~feature_importance_df.columns.str.contains(r'^Unnamed')].copy()
df

#break Category column into parts based on but keeping full info
parts = df['Category'].str.split(':', expand=True)
df['aav'] = parts[0]
df['replicate'] = parts[1].str.extract(r'replicate_(\d+)', expand=False)
df['cell_category'] = parts[2].fillna('None')
df = df.dropna(subset=['replicate'])
df['replicate'] = df['replicate'].astype(int)
df = df[df['cell_category'] != 'None']
df= df[df['aav'] != 'infection_count']

#list of all unique replicate features
rep1_features = set(
	df[(df['replicate'] == 1) & (df['cell_category']=='Excitatory') & (df['aav']=='retro')]["Feature"].unique().tolist()
)

rep2_features = set(
	df[(df['replicate'] == 2) & (df['cell_category']=='Excitatory') & (df['aav']=='retro')]["Feature"].unique().tolist()
)

rep3_features = set(
	df[(df['replicate'] == 3) & (df['cell_category']=='Excitatory') & (df['aav']=='retro')]["Feature"].unique().tolist()
)

#venn diagram
from matplotlib_venn import venn3
plt.figure(figsize=(6,6))
venn3(
	[rep1_features, rep2_features, rep3_features],
	set_labels = ('Replicate 1', 'Replicate 2', 'Replicate 3')
)

#show the features that exist in at least 2 replicates
common_features = (
	rep1_features.intersection(rep2_features)
	.union(rep1_features.intersection(rep3_features))
	.union(rep2_features.intersection(rep3_features))
)

#calculate the avg feature importance for each common feature then sort descending
common_feature_importance = {}
for feature in common_features:
	avg_importance = df[df['Feature'] == feature]['Importance'].mean()
	common_feature_importance[feature] = avg_importance

#sort the common features by avg importance
sorted_common_features = sorted(
	common_feature_importance.items(),
	key=lambda x: x[1],
	reverse=True
)
#print the sorted common features with their avg importance
for feature, importance in sorted_common_features:
	print(f"Feature: {feature}, Avg Importance: {importance:.4f}")




###########################################################################
# Boxplot of all features
###########################################################################

#open te feature importance csv
feature_importance_df = pd.read_csv("runs/feature_importances_replicates.csv")

#clean up columns
df = feature_importance_df.copy()
df = df.loc[:, ~df.columns.str.contains(r'^Unnamed')]
key_col = df.select_dtypes(include='object').columns[0]
df = df[[key_col]].rename(columns={key_col:'key'})

#parse key structure
parts = df['key'].str.split(':', expand=True)
df['aav'] = parts[0]
df['replicate'] = parts[1].str.extract(r'replicate_(\d+)', expand=False)
df['category'] = parts[2].fillna('None')
df = df.dropna(subset=['replicate'])
df['replicate'] = df['replicate'].astype(int)

df = df[df['category'] != 'None']

#count entries per (aav, replicate, category)
counts = (
    df.groupby(['aav','replicate','category'], dropna=False)
      .size()
      .reset_index(name='n')
)

#sum counts across categories per replicate
replicate_counts = (
    counts.groupby(['aav','replicate'], as_index=False)['n']
    .sum()
)

#explicitly order replicates 1→3
replicate_counts['replicate'] = pd.Categorical(replicate_counts['replicate'], categories=[1,2,3], ordered=True)

#how many for each replicate
replicate_summary = replicate_counts.groupby('replicate')['n'].describe()
print(replicate_summary)



#plot boxplots (one per replicate)
plt.figure(figsize=(6,4))
ax = sns.boxplot(data=replicate_counts, x='replicate', y='n', order=[1,2,3], showfliers=False)
sns.stripplot(
    data=replicate_counts,
    x='replicate',
    y='n',
    order=[1,2,3],
    color='black',
    alpha=0.8,
    size=5,
    jitter=True
)
ax.set_xlabel('replicate')
ax.set_ylabel('total count per AAV')
sns.despine()
plt.tight_layout()
plt.show()