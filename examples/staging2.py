import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

#open the PCCM file
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/archive/data/PCCMwithRaw.h5ad")

parameters = {
	"adataObject": adata,
	"name": "PCCM",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	#"excludeGenes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": ".",
	"outputName": "PCCM.html",
	"excludeMarkers" : True,
	"mapOnly" : False,
	"showRF": True,
	"removeOutliers": True,
	"runCellTypist": False,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"minify": False,
}

results = t.TAP(**parameters)