import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/TAP/examples/data/Cortex.h5ad",
	"name": "Cortex: M4: Kash (rep 3); M5: 1-2 (reps 1,2)",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	#"genes": ["SHH10","Retro"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicates","Cell Type"],
	"exclude" : ["Microglia"],
	"outputPath": ".",
	"outputName": "Cortex4.html",
	"excludeMarkers" : True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"minify": False,
	"rfHyperParameterTune":True,
	"rfHyperParameterIterations": 10,
	"rfPermuteFeatureImportance": True,
	"rfPermuteRepeats": 2,
}

results = t.TAP(**parameters)
