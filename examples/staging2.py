import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

#open the PCCM file
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/archive/data/PCCMwithRaw.h5ad")

parameters = {
	"adataObject": adata,
	"name": "PCCM- Test",
	"genes": ["AAV2","Retro"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": ".",
	"outputName": "PCCM.html",
	"excludeMarkers" : True,
	"mapOnly" : False,
	"showRF": True,
	"removeOutliers": True,
	"runCellTypist": True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"minify": False,
}

results = t.TAP(**parameters)