import tap as t
import scanpy as sc
import pandas as pd
import numpy as np


parameters = {
	"filename" : "../assets/data/PCCM.h5ad",
	"name": "PCCM",
	"genes": ["AAV2","Retro"],
	"excludeGenes" : ["AAV1","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Anc80","SHH10"], 
	"categories": ["cell_type"],
	"outputPath": ".", #current directory
	"outputName": "test.html",
	"runCellTypist" : False,
	"cellTypistModel": "BrainCellData_268_CellTypes.pkl",
	"excludeMarkers" : True,
	"mapOnly" : False,
	"showRF": True,
	"removeOutliers": True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"rfType": "classifier",
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 10,
	"rfPermuteFeatureImportance": False,
	"rfPermuteRepeats": 3,
	"minify": False,
}

results = t.TAP(**parameters)