import tap as t

parameters = {
	"filename" : "../assets/data/PCCM.h5ad",
	"name": "PCCM",
	"genes": ["AAV2","Retro"],
	"categories": ["cell_type"],
	"outputPath": ".", #current directory
	"outputName": "test.html",
	"runCellTypist" : False,
	"excludeMarkers" : True,
	"mapOnly" : False,
	"showRF": True,
	"removeOutliers": True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"rfHyperParameterTune": True,
	"rfHyperParameterIterations": 2,
}

results = t.TAP(**parameters)
