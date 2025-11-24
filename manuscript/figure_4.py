import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tap as t
from datetime import datetime
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from matplotlib_venn import venn3
sns.set(style="whitegrid")
sc.settings.set_figure_params(dpi=300)

# #opens the b19 dataset
# adata = sc.read_h5ad("data/B19_noRG_preprocessed.h5ad")

# #removes any null cells
# adata = adata[~adata.obs['common_name_liger_granular'].isnull()].copy()

# #saves the adata object for later use
# adata.write("data/B19_noRG_preprocessed_cleaned.h5ad")


#fresh dataset
adata = sc.read_h5ad("data/B19_noRG_preprocessed_cleaned.h5ad")
pd.options.display.max_columns = None

#next, we're going to split the index into the replicate information
#index as an obs column so we can work with it as a string
adata.obs["cell_barcode"] = adata.obs.index

#everything up the last _ in the cell barcode is the replicate
adata.obs['replicate_1'] = adata.obs['cell_barcode'].str.rsplit('_', n=1).str[0]

#everything up to the 2nd to last _ in the cell barcode is the replicate2
adata.obs['replicate_2'] = adata.obs['replicate_1'].str.rsplit('_', n=1).str[0]

#counts of each replicate "level"
adata.obs['replicate_1'].value_counts()
adata.obs['replicate_2'].value_counts()

#iterate replicate 2 and plot hitograms of total counts and viral umis
for rep in adata.obs['replicate_2'].unique():
	subset = adata[adata.obs['replicate_2'] == rep]
	fig, axs = plt.subplots(1, 2, figsize=(13, 4))
	sns.histplot(subset.obs['total_counts'], bins=50, ax=axs[0], color='blue')
	axs[0].set_title(f'Replicate {rep}\nTotal Counts')
	sns.histplot(subset.obs['viral_umis'], bins=50, ax=axs[1], color='orange')
	axs[1].set_title(f'Replicate {rep}\nViral UMIs')
	plt.tight_layout()
	plt.show()

#makes this a bit cleaner to read for the figure
adata.obs["Total Viral UMIs"] = adata.obs["viral_umis"]


#after exploring, we use the followin
replicates = ["CCEx16_B19_noRG_1x","CCEx15_B19_noRG_1x","CCEx12_B19_noRG_1in10x"]

#filter to only these replicates
adata = adata[adata.obs['replicate_2'].isin(replicates)].copy()

#add replicates column for TAP
adata.obs["replicates"] = adata.obs['replicate_2']

#rename replicates to replicate_n for simplicity
adata.obs['replicates'] = adata.obs['replicate_2'].replace({
	"CCEx16_B19_noRG_1x": "Replicate_1",
	"CCEx15_B19_noRG_1x": "Replicate_2",
	"CCEx12_B19_noRG_1in10x": "Replicate_3",
})

#save the dataset for later use
#adata.write("data/B19_noRG_preprocessed_replicates.h5ad")

#generate a TAP for replicates
parameters = {
	"filename": "data/B19_noRG_preprocessed_replicates.h5ad",
	"name": "B19 : Regressor : Replicates 1-3",
	"genes": ["N", "P", "M", "EGFP", "L","viral_umis"],
	#"genes": ["EGFP","viral_umis"],
	"excludeGenes": ["AAV_INTERGENIC"],
	"exclude": ["Oligodendrocyte","Astrocyte_Glutamatergic_Neuron","OPC","NPC","Mitotic"],
	"categories": ["replicates","common_name_liger_granular"],
	"categoryNames": ["Replicate","Cell Type"],
	"outputPath": "runs/",
	"outputName": "b19_regressor_replicates2.html",
	"clusterMethod": "gaussian", #gaussian or threshold
	"rfType": "regressor", #classifier or regressor
	"useAllGenes": False,
	"removeOutliers": False,
	"replicate_mode": True,
	"excludeMarkers" : False,
}

#build the TAP
results = t.TAP(**parameters)


#open fresh dataset
adata = sc.read_h5ad("data/B19_noRG_preprocessed_replicates.h5ad")

#remove Oligodendrocytes and Astrocyte_Glutamatergic_Neuron
adata = adata[~adata.obs['common_name_liger_granular'].isin(["Oligodendrocyte","Astrocyte_Glutamatergic_Neuron"])].copy()

#define a custom color palette
custom_colors = [
    '#0000ff',   #Blue
    '#ff6600',   #Orange
    '#008000',   #Green
    '#ff0000',   #Red
    '#b645f3',   #purple
	'#8c564b',   #brown
	'#e377c2',   #pink
]

#umap of cell types
sc.pl.umap(
    adata,
    color=["common_name_liger_granular"],
    frameon=False,
    wspace=0.4,
    hspace=0.4,
    ncols=3,
    show=True,
    title="Cell Types",
    palette=custom_colors
)

#log1p the viral umis for better visualization
adata.obs["Total Viral UMIs Log1p"] = np.log1p(adata.obs["Total Viral UMIs"])

#same for N,M,P,EGFP,L
for gene in ["N", "P", "M", "EGFP", "L"]:
	adata.obs[gene + " Log1p"] = np.log1p(adata.obs[gene])


tap_colors = [(0.0, 'blue'), (0.1, 'white'), (1.0, 'red')]
tap_cmap = LinearSegmentedColormap.from_list('tap_whitebluered', tap_colors, N=256)


#umap of cell types
sc.pl.umap(
    adata,
    color=["Total Viral UMIs Log1p"],
    frameon=False,
    wspace=0.4,
    hspace=0.4,
    ncols=3,
    show=True,
	cmap=tap_cmap,
)

#umap of N, P, M, EGFP, L
for protein in ["N", "P", "M", "EGFP", "L"]:
	sc.pl.umap(
		adata,
		color=[f"{protein} Log1p"],
		frameon=False,
		wspace=0.4,
		hspace=0.4,
		ncols=3,
		show=True,
	)



#dotplot with grouped cell types
sc.pl.dotplot(
    adata,
    var_names=["N", "P", "M", "EGFP", "L","Total Viral UMIs"],
    groupby="common_name_liger_granular",
    use_raw=False,
    dendrogram=False,
	swap_axes=True,
    show=True
)

########################################################################################
#Supplementary Figure: Feature importance counts per replicate
########################################################################################

#open te feature importance csv
feature_importance_df = pd.read_csv("runs/feature_importances_sadb19.csv")

#clean up columns
df = feature_importance_df.copy()
df = df.loc[:, ~df.columns.str.contains(r'^Unnamed')]
key_col = df.select_dtypes(include='object').columns[0]
df = df[[key_col]].rename(columns={key_col:'key'})

#parse key structure
parts = df['key'].str.split(':', expand=True)
df['aav'] = parts[0]
df['replicate'] = parts[1].str.extract(r'Replicate_(\d+)', expand=False)
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
plt.figure(figsize=(6,5))
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
ax.set_ylabel('total count per sadb19 gene')
sns.despine()
plt.tight_layout()
plt.show()


########################################################################################
#Supplementary Figure:  venn diagram of common features across replicates (GLut)
########################################################################################


#imports
import pandas as pd

#io
feature_importance_df=pd.read_csv("runs/feature_importances_sadb19.csv")
feature_importance_df

#drop stray index cols
df=feature_importance_df.loc[:, ~feature_importance_df.columns.str.contains(r'^Unnamed')].copy()
df

#break Category column into parts based on but keeping full info
parts = df['Category'].str.split(':', expand=True)
df['aav'] = parts[0]
df['replicate'] = parts[1].str.extract(r'Replicate_(\d+)', expand=False)
df['cell_category'] = parts[2].fillna('None')
df = df.dropna(subset=['replicate'])
df['replicate'] = df['replicate'].astype(int)
df = df[df['cell_category'] != 'None']
df= df[df['aav'] != 'infection_count']
df


#list of all unique replicate features
rep1_features = set(
	df[(df['replicate'] == 1) & (df['cell_category']=='Glutamatergic_Neuron') & (df['aav']=='viral_umis')]["Feature"].unique().tolist()
)

rep2_features = set(
	df[(df['replicate'] == 2) & (df['cell_category']=='Glutamatergic_Neuron') & (df['aav']=='viral_umis')]["Feature"].unique().tolist()
)

rep3_features = set(
	df[(df['replicate'] == 3) & (df['cell_category']=='Glutamatergic_Neuron') & (df['aav']=='viral_umis')]["Feature"].unique().tolist()
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



########################################################################################
#Supplementary Figure:  venn diagram of common features across replicates (Interneurons)
########################################################################################


#imports
import pandas as pd

#io
feature_importance_df=pd.read_csv("runs/feature_importances_sadb19.csv")
feature_importance_df

#drop stray index cols
df=feature_importance_df.loc[:, ~feature_importance_df.columns.str.contains(r'^Unnamed')].copy()
df

#break Category column into parts based on but keeping full info
parts = df['Category'].str.split(':', expand=True)
df['aav'] = parts[0]
df['replicate'] = parts[1].str.extract(r'Replicate_(\d+)', expand=False)
df['cell_category'] = parts[2].fillna('None')
df = df.dropna(subset=['replicate'])
df['replicate'] = df['replicate'].astype(int)
df = df[df['cell_category'] != 'None']
df= df[df['aav'] != 'infection_count']
df


#list of all unique replicate features
rep1_features = set(
	df[(df['replicate'] == 1) & (df['cell_category']=='Interneuron') & (df['aav']=='viral_umis')]["Feature"].unique().tolist()
)

rep2_features = set(
	df[(df['replicate'] == 2) & (df['cell_category']=='Interneuron') & (df['aav']=='viral_umis')]["Feature"].unique().tolist()
)

rep3_features = set(
	df[(df['replicate'] == 3) & (df['cell_category']=='Interneuron') & (df['aav']=='viral_umis')]["Feature"].unique().tolist()
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


