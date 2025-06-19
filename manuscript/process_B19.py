import scanpy as sc
import pandas as pd
import numpy as np
import tap as t
import matplotlib.pyplot as plt
import tap as t
pd.set_option('display.max_columns', None)

#find any genes that contain CVS by grepping case insensitive
viral_genes_b19 = ["N", "P", "M", "EGFP", "L"] #b19

#open the data
adata = sc.read_h5ad("data/B19_noRG.h5ad")
adata.shape

#open csv file
meta_data = pd.read_csv("data/B19_noRG_Liger_Summary_Complete_counts.csv", index_col=0)
meta_data.shape

#make the CBC column the index
meta_data.set_index("CBC", inplace=True)

#add metadata to adata matching adata index	to meta_data index
adata.obs = adata.obs.join(meta_data, how="left")
adata.obs


#add metadata
adata.var["mt"] = adata.var_names.str.startswith(("mt-","MT-","Mt-"))
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

#qc metrics using scanpy
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

#plotqc metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
adata.obs.total_counts.hist(bins=100)
plt.xticks(np.arange(0, 50000, 1000))
plt.xticks(rotation=45)

#filter cells and genes
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, max_counts=50000)
sc.pp.filter_genes(adata, min_cells=3)

#doublet detection using scrublet
sc.pp.scrublet(adata, threshold=0.5)

#save viral_umis to adata.obs
adata.obs["viral_umis"] = adata[:, viral_genes_b19].X.sum(axis=1)

#add all viral genes as observation metadata for easier access later
for gene in viral_genes_b19:
	adata.obs[gene] = adata[:, gene].X.sum(axis=1)

#save raw data
adata.raw = adata.copy()

#normalize and log
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

#save the file as
#adata.write_h5ad("data/B19_noRG_preprocessed.h5ad")

#open the preprocessed data
adata = sc.read_h5ad("data/B19_noRG_preprocessed.h5ad")
adata.shape
adata.obs

#take a look at the viral genes
adata[:,viral_genes_b19].to_df()
sc.pl.scatter(adata, "total_counts", "viral_umis", title="Viral UMIs vs Total Counts")
adata.obs["viral_umis"].hist(bins=100)

#calculate % of viral UMIs in total counts
adata.obs["viral_umis_pct"] = adata.obs["viral_umis"] / adata.obs["total_counts"]

#count how many cells have viral UMIs that are more than 10% of counts?
viral_umis_pct_threshold = 0.08
adata.obs["viral_umis_pct"].hist(bins=100)
print(f"Number of cells with viral UMIs > {viral_umis_pct_threshold*100}% of total counts: "
	  f"{(adata.obs['viral_umis_pct'] > viral_umis_pct_threshold).sum()}")

#remove cells with viral UMIs that are more than 8% of total counts
adata = adata[adata.obs["viral_umis_pct"] <= viral_umis_pct_threshold].copy()

#take a look at the viral genes again
adata[:,viral_genes_b19].to_df()
sc.pl.scatter(adata, "total_counts", "viral_umis", title="Viral UMIs vs Total Counts")
adata.obs["viral_umis"].hist(bins=100)

#save the file as
#adata.write_h5ad("data/B19_noRG_preprocessed.h5ad")

#open the preprocessed data
adata = sc.read_h5ad("data/B19_noRG_preprocessed.h5ad")
adata.shape
adata.obs

#remove the viral genes from the adata object in a memory efficient way
adata = adata[:, ~adata.var_names.isin(viral_genes_b19)].copy()

#highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pl.highly_variable_genes(adata)

#get list of highly variable genes
hvg_genes = adata.var_names[adata.var["highly_variable"]].tolist()

#create a filter list with highly variable genes and viral genes
filter_list = hvg_genes

#remove duplicates from the filter list
filter_list = list(set(filter_list))

#filter the adata
adata.var_names_make_unique()
adata = adata[:, filter_list]
adata.shape

#pca
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

#clustering & umap
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=0.2)
sc.tl.umap(adata)

#take a look at the metadata
adata.obs

#cast the ledien column to string
adata.obs["leiden"] = adata.obs["leiden"].astype(str)
adata.obs["leiden"] = "leiden_"+adata.obs["leiden"]
adata.obs["leiden"].value_counts()

#umap plots
sc.pl.umap(adata)
sc.pl.umap(adata,color=["viral_umis"])
sc.pl.umap(adata,color=["common_name_liger_coarse"])

#save the adata object
#adata.write_h5ad("data/B19_noRG_preprocessed.h5ad")

#open the preprocessed data
adata = sc.read_h5ad("data/B19_noRG_preprocessed.h5ad")

marker_genes = ["Gad1", "Gad2", "Slc32a1", "Slc17a7", "Slc17a6", "Sncg", "Olig2",
				"Vip", "Sst", "Chodl", "Pvalb", "Rorb", "Fezf2", "Sulf1",
				"Foxp2", "Nxph4", "Aqp4", "Mbp", "Cldn5", "Ctss", "C1qa"]

#remove any genes that are not in the adata object
marker_genes = [gene for gene in marker_genes if gene in adata.var_names]

#explore the obs before running TAP
adata = sc.read_h5ad("data/B19_noRG_preprocessed.h5ad")
adata.obs


#run for classifier
parameters = {
	"filename": "data/B19_noRG_preprocessed.h5ad",
	"name": "B19 : Classifier",
	"genes": ["N", "P", "M", "EGFP", "L","viral_umis"],
	"excludeGenes": ["AAV_INTERGENIC"],
	"exclude": ["NA","Oligodendrocyte","Astrocyte_Glutamatergic_Neuron",""],
	"useAllGenes": False,
	"categories": ["common_name_liger_granular"],
	"categoryNames": ["Cell Type"],
	"outputPath": ".",
	"outputName": "b19_classifier.html",
	"excludeMarkers" : False,
	"mapOnly" : False,
	"clusterMethod": "gaussian",
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 10,
	"rfType": "classifier", #classifier or regressor
	"balance": "smote", #"smote" or None
	"minify": False,
	"removeOutliers": False,
	"replicate_mode": False,
}

results = t.TAP(**parameters)


#run for regressor
parameters = {
	"filename": "data/B19_noRG_preprocessed.h5ad",
	"name": "B19 : Regressor",
	"genes": ["N", "P", "M", "EGFP", "L","viral_umis"],
	"excludeGenes": ["AAV_INTERGENIC"],
	"exclude": ["NA","Oligodendrocyte","Astrocyte_Glutamatergic_Neuron",""],
	"useAllGenes": False,
	"categories": ["common_name_liger_granular"],
	"categoryNames": ["Cell Type"],
	"outputPath": ".",
	"outputName": "b19_classifier.html",
	"excludeMarkers" : False,
	"mapOnly" : False,
	"clusterMethod": "gaussian",
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 10,
	"rfType": "regressor", #classifier or regressor
	"balance": "smote", #"smote" or None
	"minify": False,
	"removeOutliers": False,
	"replicate_mode": False,
}

results = t.TAP(**parameters)
