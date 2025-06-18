import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from biothings_client import get_client

#open the h5ad file
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj.Robj.h5ad")
adata.obs["replicate"] = np.random.randint(0, 2, size=adata.obs.shape[0])
adata.obs["replicate"] = adata.obs["replicate"].astype(str)
adata.obs["replicate"] = "replicate_" + adata.obs["replicate"]
adata.obs.groupby("replicate")["xincelltype230416"].value_counts()

parameters = {
	"adataObject": adata,
	"name": "Xin Jin 2024 - Test 2",
	"genes": ["BC1a","BC1b","BC1c",],
	"exclude" : ["Microglia","Mural","Fibroblast"],
	"categories": ["replicate","xincelltype230416"],
	"categoryNames": ["Replicate","Cell Type"],
	"outputPath": ".",
	"outputName": "xin_test2.html",
	"excludeMarkers" : True,
	"mapOnly" : False,
	"showRF": True,
	"removeOutliers": True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"rfType": "classifier",
	"balance": "smote",
	"minify": False,
}

results = t.TAP(**parameters)