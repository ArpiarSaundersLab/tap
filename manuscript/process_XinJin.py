import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from biothings_client import get_client
import seaborn as sns
from datetime import datetime

#open the h5ad file
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj_2.Robj.h5ad")
# adata.obs["replicate"] = np.random.randint(0, 2, size=adata.obs.shape[0])
# adata.obs["replicate"] = adata.obs["replicate"].astype(str)
# adata.obs["replicate"] = "replicate_" + adata.obs["replicate"]
# adata.obs.groupby("replicate")["xincelltype230416"].value_counts()


parameters = {
	"adataObject": adata,
	"name": "Xin Jin 2024 | run:"+datetime.now().strftime("%d-%m-%Y"),
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
	"excludeGenes": [
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
	#"genes": ["BC1a","BC1b","BC1c","BC2a","BC2b","BC2c","BC4a","BC4b","BC4c","BC10a","BC10b","BC10c"],
	"genes": ["BC1a","BC1b","BC1c","BC2a","BC2b","BC2c","BC4a","BC4b","BC4c"],
	"exclude" : ["Microglia","Mural","Fibroblast","CajalRetzius cells","IN nonMGE"],
	"categories": ["xincelltype230416"],
	"categoryNames": ["Cell Type"],
	"outputPath": "runs/",
	"outputName": "xin_jin_"+str(datetime.now().strftime("%d-%m-%Y"))+".html",
	"mapOnly" : False,
	"showRF": True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"rfType": "classifier",
	"runCellTypist": False,
	"cellTypistModel": "Developing_Mouse_Brain.pkl",
	"cellTypistLevels": 1,
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 20,
	"rfPermuteFeatureImportance": False,
	"rfPermuteRepeats": 5,
	"balance": "smote", #None or "smote"
	"minify": False,
	"removeOutliers": True,
	"replicate_mode": True,
	"excludeMarkers" : True,
	"useAllGenes": False,
	"minCells":20, 
	"minSeroTypeCells":20,

}

results = t.TAP(**parameters)