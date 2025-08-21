import tap as t

#Cortex M6
parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/molotAAV/Code/Analyses/M6/Mouse/Cortex/Data/Dialout/preprocessed.h5ad",
	"name": "Cortex: M6",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": "runs/",
	"outputName": "Cortex_M6.html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"replicate_mode": False,
	"excludeMarkers" : True,
}
results = t.TAP(**parameters)


#Midbrain M6
parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/molotAAV/Code/Analyses/M6/Mouse/Midbrain/Data/Dialout/preprocessed.h5ad",
	"name": "Midbrain: M6",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"outputPath": "runs/",
	"outputName": "Midbrain_M6.html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"replicate_mode": False,
	"excludeMarkers" : True,
}
results = t.TAP(**parameters)