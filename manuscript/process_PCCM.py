import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

#open the PCCM file
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/archive/data/PCCMwithRaw.h5ad")

#add C1ql3 + Retro, save as obs c1_retro
#idea is to try and find a o receptor like this.
genes = ["Retro", "C1ql3","Camk2d"]
gene_indices = [adata.raw.var_names.get_loc(g) for g in genes]

X = adata.raw.X[:, gene_indices]
if not isinstance(X, np.ndarray):
    X = X.toarray()

# mask where both expressions are > 0
mask = (X[:, 0] > 0) & (X[:, 1] > 0)

# sum only where both > 0, else 0
c1_retro = np.where(mask, X.sum(axis=1), 0)

adata.obs["c1_retro"] = c1_retro
adata[:,"Retro"].to_df()
adata[:,"C1ql3"].to_df()
adata.obs["c1_retro"].hist(bins=100)




parameters = {
	"adataObject": adata,
	"name": "PCCM",
	#"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"genes": ["AAV2","c1_retro"],
	"excludeGenes": ["Camk2d","C1ql3","AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": ".",
	"outputName": "data/PCCM2.html",
	"excludeMarkers" : True,
	"removeOutliers": False,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 6,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"minify": False,
}

results = t.TAP(**parameters)