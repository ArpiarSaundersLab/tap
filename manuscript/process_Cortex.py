import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad",
	"name": "Cortex: M4: Kash (rep 3); M5: 1-2 (reps 1,2)",
	"genes": ["AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80"],
	"excludeGenes": ["Retro_Kash","Retro_10","SHH10","AAV1"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	"exclude" : ["Microglia","NPCs","OPCs"],
	"outputPath": "runs/",
	"outputName": "Cortex.html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 5,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"replicate_mode": True,
	"excludeMarkers" : True,
}

results = t.TAP(**parameters)
