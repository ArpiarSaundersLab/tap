import tap as t

#Cortex M3
parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/molotAAV/Code/Analyses/M3/Data/preprocessed.h5ad",
	"name": "Cortex: M3",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": "runs/",
	"outputName": "Cortex_M3.html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"replicate_mode": False,
	"excludeMarkers" : True,
}
results = t.TAP(**parameters)