import scanpy as sc
import pandas as pd
import numpy as np
import tap as t
import matplotlib.pyplot as plt
import tap as t
import anndata as ad
import scanpy.external as sce
pd.set_option('display.max_columns', None)

#find any genes that contain CVS by grepping case insensitive
viral_genes_cvs  = ["N", "P", "M", "tdTom", "L"] #cvs

#open the data and store and object to list
adata_object = []

#iterate all files and add to the anndata object
for i in range(1,6):
	print(f"imported file {i}")

	#open the data
	adata = sc.read_h5ad(f"data/CCEx17_CVS_noRG_1x_L{i}.h5ad")
	
	#add metadata
	adata.var["mt"] = adata.var_names.str.startswith(("mt-","MT-","Mt-"))
	adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
	adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

	#qc metrics using scanpy
	sc.pp.calculate_qc_metrics(
		adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
	)

	#plotqc metrics
	adata.obs.total_counts.hist(bins=100)

	#filter cells and genes
	sc.pp.filter_cells(adata, min_counts=2000)
	sc.pp.filter_cells(adata, max_counts=35000)
	sc.pp.filter_genes(adata, min_cells=3)	

	#plotqc metrics
	adata.obs.total_counts.hist(bins=100)

	adata_object.append(adata)
	adata_object[i-1].obs["replicate"] = f"Replicate_L{i}"

#take a look at the list of objects
adata_object

#concatenate all together
adata = ad.concat(adata_object)

#look at the full object
adata
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
adata.obs.total_counts.hist(bins=100)

#save viral_umis to adata.obs
adata.obs["viral_umis"] = adata[:, viral_genes_cvs].X.sum(axis=1)

#add all viral genes as observation metadata for easier access later
for gene in viral_genes_cvs:
	adata.obs[gene] = adata[:, gene].X.sum(axis=1)

#save raw data
adata.raw = adata.copy()

#normalize and log
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

#save the file as
#adata.write_h5ad("data/CCEx17_CVS_noRG_all_preprocessed.h5ad")

#open the preprocessed data
adata = sc.read_h5ad("data/CCEx17_CVS_noRG_all_preprocessed.h5ad")
adata.shape
adata.obs

#take a look at the viral genes
adata[:,viral_genes_cvs].to_df()
sc.pl.scatter(adata, "total_counts", "viral_umis", title="Viral UMIs vs Total Counts")
adata.obs["viral_umis"].hist(bins=100)

#calculate % of viral UMIs in total counts
adata.obs["viral_umis_pct"] = adata.obs["viral_umis"] / adata.obs["total_counts"]

#count how many cells have viral UMIs that are more than 10% of counts?
viral_umis_pct_threshold = 0.10
adata.obs["viral_umis_pct"].hist(bins=100)
print(f"Number of cells with viral UMIs > {viral_umis_pct_threshold*100}% of total counts: "
	  f"{(adata.obs['viral_umis_pct'] > viral_umis_pct_threshold).sum()}")

#remove cells with viral UMIs that are more than 8% of total counts
adata = adata[adata.obs["viral_umis_pct"] <= viral_umis_pct_threshold].copy()

#take a look at the viral genes again
adata[:,viral_genes_cvs].to_df()
sc.pl.scatter(adata, "total_counts", "viral_umis", title="Viral UMIs vs Total Counts")
adata.obs["viral_umis"].hist(bins=100)

#save the file as
#adata.write_h5ad("data/CCEx17_CVS_noRG_all_preprocessed.h5ad")

#open the preprocessed data
adata = sc.read_h5ad("data/CCEx17_CVS_noRG_all_preprocessed.h5ad")
adata.shape
adata.obs

#remove the viral genes from the adata object in a memory efficient way
adata = adata[:, ~adata.var_names.isin(viral_genes_cvs)].copy()

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
adata.obs
#pca
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

#harmony integration
sce.pp.harmony_integrate(adata, key='replicate')


#clustering & umap
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25, use_rep='X_pca_harmony')
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

#save the adata object
#adata.write_h5ad("data/CCEx17_CVS_noRG_all_preprocessed.h5ad")

#open the preprocessed data
adata = sc.read_h5ad("data/CCEx17_CVS_noRG_all_preprocessed.h5ad")

marker_genes = ["Gad1", "Gad2", "Slc32a1", "Slc17a7", "Slc17a6", "Sncg", "Olig2",
				"Vip", "Sst", "Chodl", "Pvalb", "Rorb", "Fezf2", "Sulf1",
				"Foxp2", "Nxph4", "Aqp4", "Mbp", "Cldn5", "Ctss", "C1qa"]

#remove any genes that are not in the adata object
marker_genes = [gene for gene in marker_genes if gene in adata.var_names]


#plot marker genes
sc.pl.umap(adata, color=marker_genes, ncols=5, wspace=0.4, hspace=0.4,
		   frameon=False, title=marker_genes, use_raw=False)
sc.pl.umap(adata, color="leiden", frameon=False, title="leiden clusters")

cluster2annotations = {
	"leiden_0": "Glutamatergic neurons",
	"leiden_1": "OPCs",
	"leiden_2": "OPCs",
	"leiden_3": "Gabaergic neurons",
	"leiden_4": "Glutamatergic neurons",
	"leiden_5": "Astrocytes",
	"leiden_6": "Oligodendrocytes",
	"leiden_7": "Oligodendrocytes",
	"leiden_8": "Oligodendrocytes",
	"leiden_9": "Gabaergic neurons",
	"leiden_10": "Microglia"
}


#add the annotations to the adata object
adata.obs["cell_type"] = adata.obs["leiden"].map(cluster2annotations)

#plot the cell types
sc.pl.umap(adata, color="cell_type")



#save the adata object
#adata.write_h5ad("data/CCEx17_CVS_noRG_all_preprocessed.h5ad")


#explore the obs before running TAP
adata = sc.read_h5ad("data/CCEx17_CVS_noRG_all_preprocessed.h5ad")
adata.obs


#run tap 
import tap as t

parameters = {
	"filename": "data/CCEx17_CVS_noRG_all_preprocessed.h5ad",
	"name": "CVS : Regressor",
	"genes": ["N", "P", "M", "tdTom", "L","viral_umis"],
	"excludeGenes": ["AAV_INTERGENIC"],
	"exclude": ["NA","Oligodendrocyte","Astrocyte_Glutamatergic_Neuron","Microglia",""],
	"useAllGenes": False,
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": ".",
	"outputName": "cvs_regressor.html",
	"excludeMarkers" : False,
	"mapOnly" : False,
	"clusterMethod": "gaussian",
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 10,
	"rfType": "regressor", #classifier or regressor
	"balance":None, #"smote" or None
	"minify": False,
	"removeOutliers": False,
	"replicate_mode": False,
}

results = t.TAP(**parameters)



parameters = {
	"filename": "data/CCEx17_CVS_noRG_all_preprocessed.h5ad",
	"name": "CVS : Classifier",
	"genes": ["N", "P", "M", "tdTom", "L","viral_umis"],
	"excludeGenes": ["AAV_INTERGENIC"],
	"exclude": ["NA","Oligodendrocyte","Astrocyte_Glutamatergic_Neuron","Microglia",""],
	"useAllGenes": False,
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": ".",
	"outputName": "cvs_cegressor.html",
	"excludeMarkers" : False,
	"mapOnly" : False,
	"clusterMethod": "gaussian",
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 10,
	"rfType": "regressor", #classifier or regressor
	"balance":"smote", #"smote" or None
	"minify": False,
	"removeOutliers": False,
	"replicate_mode": False,
}

results = t.TAP(**parameters)
