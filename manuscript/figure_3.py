import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tap as t
from datetime import datetime
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
sns.set(style="whitegrid")
sc.settings.set_figure_params(dpi=300)

#opens the xinjin dataset
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj_2.Robj.h5ad")

#generates a TAP of the AAV screen in the Zheng et al 2024 paper
parameters = {
	"adataObject": adata,
	"name": "Xin Jin 2024 | run:"+datetime.now().strftime("%d-%m-%Y"),
	"genes": ["BC1a","BC1b","BC1c","BC2a","BC2b","BC2c","BC4a","BC4b","BC4c"], 
	"excludeGenes": [ #remove other AAVs from introducing bias to the results
			  "BC3a","BC3b","BC3c",
			  "BC5a","BC5b","BC5c",
			  "BC6a","BC6b","BC6c",
			  "BC7a","BC7b","BC7c",
			  "BC8a","BC8b","BC8c",
			  "BC9a","BC9b","BC9c",
			  "BC11a","BC11b","BC11c",
			  "BC12a","BC12b","BC12c",
			  "BC13a","BC13b","BC13c",
			  "BC14a","BC14b","BC14c"],
	"exclude" : ["Microglia","Mural","Fibroblast","CajalRetzius cells","IN nonMGE"],
	"categories": ["xincelltype230416"],
	"categoryNames": ["Cell Type"],
	"outputPath": "runs/",
	"outputName": "xin_jin_"+str(datetime.now().strftime("%d-%m-%Y"))+".html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"rfType": "classifier",
	"balance": "smote",
	"minify": False,
	"replicate_mode": True,
	"excludeMarkers" : True,
	"useAllGenes": False,
	"minCells":20, 
	"minSeroTypeCells":20
}

results = t.TAP(**parameters)


#cell type umap with custom colors
sc.pl.umap(
    adata,
    color=["xincelltype230416"],
    frameon=False,
    wspace=0.4,
    hspace=0.4,
    ncols=3,
    show=True,
	title="Cell Types"
)


#filters unwanted cells  ["Microglia","Mural","Fibroblast","CajalRetzius cells","IN nonMGE"]
adata = adata[~adata.obs["xincelltype230416"].isin(["Microglia","Mural","Fibroblast","CajalRetzius cells","IN nonMGE"])].copy()

#order the categories in reverse of their current order
adata.obs["xincelltype230416"] = pd.Categorical(adata.obs["xincelltype230416"], categories=adata.obs["xincelltype230416"].cat.categories[::-1], ordered=True)

#plots with grouped cell types
sc.pl.dotplot(
    adata,
    var_names=["BC1a","BC1b","BC1c","BC2a","BC2b","BC2c","BC4a","BC4b","BC4c"],
    groupby="xincelltype230416",
    use_raw=False,
    dendrogram=False,
    show=True
)

#loads fresh data
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj_2.Robj.h5ad")

#add BC1a + BC1b + BC1c expression
adata.obs["BC1_total"] = adata.obs[["BC1a", "BC1b", "BC1c"]].sum(axis=1)

#converts to log1p scale
adata.obs["BC1_total"] = np.log1p(adata.obs["BC1_total"])


tap_colors = [(0.0, 'blue'), (0.1, 'white'), (1.0, 'red')]
tap_cmap = LinearSegmentedColormap.from_list('tap_whitebluered', tap_colors, N=256)

#umap of BC1_total expression
sc.pl.umap(
	adata,
	color=["BC1_total"],
	frameon=False,
	title="BC1 (a+b+c) Expression",
	wspace=0.4,
	hspace=0.4,
	ncols=2,
	show=True,
	use_raw=False,
	s=10,
	cmap=tap_cmap
)

#umap of Lrfn5 expression
sc.pl.umap(
	adata,
	color=["Lrfn5"],
	frameon=False,
	title="Lrfn5 Expression",
	wspace=0.4,
	hspace=0.4,
	ncols=2,
	show=True,
	use_raw=False,
	s=10
)

#umap of Hs3st4 expression
sc.pl.umap(
	adata,
	color=["Hs3st4"],
	frameon=False,
	title="Lrfn5 Expression",
	wspace=0.4,
	hspace=0.4,
	ncols=2,
	show=True,
	use_raw=False,
	s=10
)


#volin plots of bc1a, bc2a, and bc4a expression
genes_of_interest = ["BC1a", "BC1b", "BC1c"]

#loads fresh data
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj_2.Robj.h5ad")

#filters to ULPNs only
adata = adata[adata.obs["xincelltype230416"].isin(["ULPN"])].copy()

#simple outlier removal for better visualization
#the genes exist in the obs layer
for gene in genes_of_interest:
	upper_bound = np.percentile(adata.obs[gene], 99)
	adata = adata[adata.obs[gene] <= upper_bound].copy()

#volin plots - all in one plot with same scale
sc.settings.set_figure_params(dpi=300, 
figsize=(5, 6), 
fontsize=25) 

sc.pl.violin(
	adata,
	keys="BC1c",
	groupby="xincelltype230416",
	stripplot=False,
	show=True,
	use_raw=True,
	palette="gray"
)

#loads fresh data
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj_2.Robj.h5ad")

#select only ULPNs
adata = adata[adata.obs["xincelltype230416"].isin(["ULPN"])].copy()


#counts per barcode and stacked bars for BC1a/b/c
bc1_list = [c for c in ["BC1a","BC1b","BC1c"] if c in adata.obs.columns]
n_total = adata.n_obs

gene= "Lrfn5"
gene_thresh = 0.5
bc1_thresh = 0.5

#precompute gene positivity once
is_gene = (adata[:, gene].X.toarray().ravel() if hasattr(adata[:, gene].X, "toarray")
           else np.asarray(adata[:, gene].X).ravel()) > gene_thresh

neg_counts, bc_only_counts, both_counts = [], [], []
for bc in bc1_list:
    is_bc = adata.obs[bc].values > bc1_thresh
    n_bc = int(is_bc.sum())
    n_both = int((is_bc & is_gene).sum())
    neg_counts.append(n_total - n_bc)          #ulpn & bc1−
    bc_only_counts.append(n_bc - n_both)       #ulpn & bc1+ & lrfn5−
    both_counts.append(n_both)                 #ulpn & bc1+ & lrfn5+

#plot
spacing = 0.7                     # < 1.0 brings bars closer
x = np.arange(len(bc1_list)) * spacing
width = spacing * 0.9             # bars fill most of the available slot
plt.margins(x=0.02) 
plt.figure(figsize=(9,5))
b1 = plt.bar(x, neg_counts, width, color="lightgray", edgecolor="black", label="BC1−")
b2 = plt.bar(x, bc_only_counts, width, bottom=neg_counts, color="#DBF1A1", edgecolor="black",
             label="BC1+ & Lrfn5−")
b3 = plt.bar(x, both_counts, width,
             bottom=np.array(neg_counts)+np.array(bc_only_counts),
             color="#F5A3AD", edgecolor="black", label="BC1+ & Lrfn5+")

#annotate percents inside segments
for i in range(len(bc1_list)):
    bottoms = 0
    for c in [neg_counts[i], bc_only_counts[i], both_counts[i]]:
        if c>0 and n_total>0:
            plt.text(x[i], bottoms + c/2, f"{c/n_total*100:.1f}%", ha="center", 
			va="center", fontsize=20, fontweight="bold")
        bottoms += c

plt.xticks(x, bc1_list)
plt.ylabel("cells")
plt.title(f"ULPN Composition")
plt.ylim(0, max(n_total+100, 1))
plt.legend(frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.show()


#loads fresh data
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj_2.Robj.h5ad")

#select only ULPNs
adata = adata[adata.obs["xincelltype230416"].isin(["Migrating neurons"])].copy()


#counts per barcode and stacked bars for BC1a/b/c
bc1_list = [c for c in ["BC1a","BC1b","BC1c"] if c in adata.obs.columns]
n_total = adata.n_obs

#precompute gene positivity once
is_gene = (adata[:, gene].X.toarray().ravel() if hasattr(adata[:, gene].X, "toarray")
           else np.asarray(adata[:, gene].X).ravel()) > gene_thresh

neg_counts, bc_only_counts, both_counts = [], [], []
for bc in bc1_list:
    is_bc = adata.obs[bc].values > bc1_thresh
    n_bc = int(is_bc.sum())
    n_both = int((is_bc & is_gene).sum())
    neg_counts.append(n_total - n_bc)          #ulpn & bc1−
    bc_only_counts.append(n_bc - n_both)       #ulpn & bc1+ & lrfn5−
    both_counts.append(n_both)                 #ulpn & bc1+ & lrfn5+

#plot
spacing = 0.7                     # < 1.0 brings bars closer
x = np.arange(len(bc1_list)) * spacing
width = spacing * 0.9             # bars fill most of the available slot
plt.margins(x=0.02) 
plt.figure(figsize=(9,5))
b1 = plt.bar(x, neg_counts, width, color="lightgray", edgecolor="black", label="BC1−")
b2 = plt.bar(x, bc_only_counts, width, bottom=neg_counts, color="#DBF1A1", edgecolor="black",
             label="BC1+ & Hs3st4-")
b3 = plt.bar(x, both_counts, width,
             bottom=np.array(neg_counts)+np.array(bc_only_counts),
             color="#F5A3AD", edgecolor="black", label="BC1+ & Hs3st4+")

#annotate percents inside segments
for i in range(len(bc1_list)):
    bottoms = 0
    for c in [neg_counts[i], bc_only_counts[i], both_counts[i]]:
        if c>0 and n_total>0:
            plt.text(x[i], bottoms + c/2, f"{c/n_total*100:.1f}%", 
			ha="center", va="center", fontsize=20, fontweight="bold")
        bottoms += c

plt.xticks(x, bc1_list)
plt.ylabel("cells")
plt.title(f"Migrating Neurons Composition")
plt.ylim(0, max(n_total+50, 1))
plt.legend(frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.show()


##############################################################################
# Supplementary Figures
##############################################################################

#loads fresh data
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj_2.Robj.h5ad")

#exvlude from analysis
exclude = ["Microglia","Mural","Fibroblast","CajalRetzius cells","IN nonMGE"]

#list of individual barcodes in analysis
genes_of_interest = ["BC1a","BC1b","BC1c","BC2a","BC2b","BC2c","BC4a","BC4b","BC4c"]

#stacked_violin for each replicate
plt.figure(figsize=(10,5))
sc.pl.stacked_violin(
	adata,
	var_names=genes_of_interest,
	groupby="xincelltype230416",
	use_raw=False,
	dendrogram=True,
	swap_axes=True,
	show=True
)

#log1p transform the genes in obs
for gene in genes_of_interest:
	adata.obs[gene] = np.log1p(adata.obs[gene])

#exclude unwanted cell types
adata = adata[~adata.obs["xincelltype230416"].isin(exclude)].copy()

#make a heatmap of the AAV expression in cell types
sc.pl.heatmap(
	adata,
	var_names=genes_of_interest,
	groupby="xincelltype230416",
	use_raw=False,
	show=True,
	swap_axes=True,
	dendrogram=False,
	figsize=(12, 6),
)




###########################################################################
# Boxplot of all features
###########################################################################

#open te feature importance csv
feature_importance_df = pd.read_csv("runs/feature_importances_replicates_xinjin.csv")
feature_importance_df

#clean up columns
df = feature_importance_df.copy()
df = df.loc[:, ~df.columns.str.contains(r'^Unnamed')]
key_col = df.select_dtypes(include='object').columns[0]
df = df[[key_col]].rename(columns={key_col:'key'})

#parse key structure
df["Feature"] = feature_importance_df["Feature"]
df["Importance"] = feature_importance_df["Importance"]
parts = df['key'].str.split(':', expand=True)
df['aav'] = parts[0]
df = df.dropna(subset=['aav'])
df = df[~df['aav'].str.contains('infection_count')]
df['replicate'] = parts[0].str[-1].map({'a':1, 'b':2, 'c':3})
df['category'] = parts[1]
df = df.dropna(subset=['category'])
df = df[df['category'] != 'None']
df['replicate'] = df['replicate'].astype(int)
df["aav"] = df["aav"].str[:-1]
df

#count entries per (aav, replicate, category)
counts = (
    df.groupby(['aav','replicate','category'], dropna=False)
      .size()
      .reset_index(name='n')
)

#boxplot of with n per (aav, replicate)
plt.figure(figsize=(8,6))
ax = sns.boxplot(data=counts, x='replicate', y='n', showfliers=False)
sns.stripplot(
	data=counts,
	x='replicate',
	y='n',
	color='black',
	alpha=0.8,
	size=5,
	jitter=True
)



###################################################################################
# venn diagram of overlapping features
###################################################################################


#open te feature importance csv
feature_importance_df = pd.read_csv("runs/feature_importances_replicates_xinjin.csv")
feature_importance_df

#clean up columns
df = feature_importance_df.copy()
df = df.loc[:, ~df.columns.str.contains(r'^Unnamed')]
key_col = df.select_dtypes(include='object').columns[0]
df = df[[key_col]].rename(columns={key_col:'key'})

#parse key structure
df["Feature"] = feature_importance_df["Feature"]
df["Importance"] = feature_importance_df["Importance"]
parts = df['key'].str.split(':', expand=True)
df['aav'] = parts[0]
df = df.dropna(subset=['aav'])
df = df[~df['aav'].str.contains('infection_count')]
df['replicate'] = parts[0].str[-1].map({'a':1, 'b':2, 'c':3})
df['category'] = parts[1]
df = df.dropna(subset=['category'])
df = df[df['category'] != 'None']
df['replicate'] = df['replicate'].astype(int)
df["aav"] = df["aav"].str[:-1]
df

##########
# ULPNS
##########

#list of all unique replicate features
rep1_features = set(
	df[(df['replicate'] == 1) & (df['category']=='ULPN') & (df['aav']=='bc1')]["Feature"].unique().tolist()
)

rep2_features = set(
	df[(df['replicate'] == 2) & (df['category']=='ULPN') & (df['aav']=='bc1')]["Feature"].unique().tolist()
)

rep3_features = set(
	df[(df['replicate'] == 3) & (df['category']=='ULPN') & (df['aav']=='bc1')]["Feature"].unique().tolist()
)

#venn diagram
from matplotlib_venn import venn3
plt.figure(figsize=(6,6))
venn3(
	[rep1_features, rep2_features, rep3_features],
	set_labels = ('Replicate 1', 'Replicate 2', 'Replicate 3')
)

#show the features that exist in all three replicates
common_features = rep1_features.intersection(rep2_features).intersection(rep3_features)
common_features

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

####################
# Migrating Neurons
####################

#list of all unique replicate features
rep1_features = set(
	df[(df['replicate'] == 1) & (df['category']=='Migrating neurons') & (df['aav']=='bc1')]["Feature"].unique().tolist()
)

rep2_features = set(
	df[(df['replicate'] == 2) & (df['category']=='Migrating neurons') & (df['aav']=='bc1')]["Feature"].unique().tolist()
)

rep3_features = set(
	df[(df['replicate'] == 3) & (df['category']=='Migrating neurons') & (df['aav']=='bc1')]["Feature"].unique().tolist()
)

#venn diagram
from matplotlib_venn import venn3
plt.figure(figsize=(6,6))
venn3(
	[rep1_features, rep2_features, rep3_features],
	set_labels = ('Replicate 1', 'Replicate 2', 'Replicate 3')
)

#show the features that exist in all three replicates
common_features = rep1_features.intersection(rep2_features).intersection(rep3_features)
common_features

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
