import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad


parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Midbrain.h5ad",
	"name": "Midbrain: M4: 1; M5: 2-3",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicates","Cell Type"],
	"exclude" : ["replicate_1"],
	"exclude" : ["replicate_1","Microglia","NPCs","OPCs"],
	"outputPath": ".",
	"outputName": "data/Midbrain.html",
	"excludeMarkers" : True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"minify": False,
	"replicate_mode": True,
}

results = t.TAP(**parameters)
