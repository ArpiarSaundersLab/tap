import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

#open the PCCM file
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/Tbr1/all_SCT.h5ad")
adata.obs

#rename nUMI to total_counts
adata.obs["total_counts"] = adata.obs["nUMI"]

#select only the genotypes HOM and WT
adata = adata[adata.obs.genotype.isin(["HOM","WT"])]
print(adata.obs.genotype.value_counts())

#does Tbr1 exist in the dataset? If so, how do the counts look?
adata[:,"Tbr1"].to_df().hist(bins=100)

#save a raw copy of the data
adata.raw = adata.copy()

#normalize and log the data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

#take the top n most variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

#take a look at the umap that was produced in the original analysis
sc.pl.umap(adata, color=["predicted_cell_type", "genotype","Tbr1"], wspace=0.4)

#downsample the dataset to 10k cells
adata = adata[adata.obs.sample(frac=10000/adata.n_obs, random_state=1).index]

#save the adata object
adata.write_h5ad("/home/kenny/Documents/OHSU/Projects/Tbr1/Tbr1_WT_HOM_20k.h5ad")

#run tap. let's see how it does
parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/Tbr1/Tbr1_WT_HOM_20k.h5ad",
	"name": "Tbr1: WT vs HOM",
	"genes": ["Tbr1"],
	"categories": ["predicted_cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": "/home/kenny/Documents/OHSU/Projects/Tbr1/",
	"outputName": "TapTbr1WtHom.html",
	"excludeMarkers" : True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 4,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler",
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 10,
}

results = t.TAP(**parameters)