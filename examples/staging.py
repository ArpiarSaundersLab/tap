import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad


parameters = {
	"filename" : "../archive/data/GSE249416/GSE249416_AAV_ctxobj.Robj.h5ad",
	#"adataObject" : adata,
	"name": "Xin Jin 2024 - All",
	# "genes": ["BC1a","BC1b","BC1c",
	# 		  "BC2a","BC2b","BC2c",
	# 		  "BC3a","BC3b","BC3c",
	# 		  "BC4a","BC4b","BC4c",
	# 		  "BC5a","BC5b","BC5c",
	# 		  "BC6a","BC6b","BC6c",
	# 		  "BC7a","BC7b","BC7c",
	# 		  "BC8a","BC8b","BC8c",
	# 		  "BC9a","BC9b","BC9c",
	# 		  "BC10a","BC10b","BC10c",
	# 		  "BC11a","BC11b","BC11c",
	# 		  "BC12a","BC12b","BC12c",
	# 		  "BC13a","BC13b","BC13c",
	# 		  "BC14a","BC14b","BC14c",],
	"genes": ["BC1a","BC1b","BC1c"],
	"exclude" : ["Microglia","Mural","Fibroblast"],
	"categories": ["xincelltype230416"],
	"categoryNames": ["Cell Type"],
	"outputPath": ".",
	"outputName": "xin_all3.html",
	"excludeMarkers" : True,
	"mapOnly" : False,
	"showRF": True,
	"removeOutliers": True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"rfType": "classifier",
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 10,
	"rfPermuteFeatureImportance": False,
	"rfPermuteRepeats": 3,
	"minify": True,
}

results = t.TAP(**parameters)

