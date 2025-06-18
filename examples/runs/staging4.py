import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/TAP/examples/data/Midbrain.h5ad",
	"name": "Midbrain: M4: 1; M5: 2-3",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	#"genes": ["AAV1","AAV2","Retro"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicates","Cell Type"],
	"exclude" : ["replicate_1"],
	"outputPath": ".",
	"outputName": "Midbrain-M4M5.html",
	"excludeMarkers" : True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"minify": False,
}

results = t.TAP(**parameters)
