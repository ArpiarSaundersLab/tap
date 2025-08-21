import tap as t

#replicates and cell types
parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad",
	"name": "Cortex Combined: M4: Kash (rep 3, sum); M5: 1-2 (reps 1,2); M6: 1",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	#"exclude" : ["Microglia","NPCs","OPCs","Astrocytes"],
	"exclude" : ["Microglia","NPCs","OPCs"],
	"outputPath": "runs/",
	"outputName": "Cortex_full.html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 5,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"replicate_mode": True,
	"excludeMarkers" : True,
}
results = t.TAP(**parameters)


#excitatory and inhibitory only
import tap as t
parameters = {
	"filename": "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad",
	"name": "Cortex Combined:  M4: Kash (rep 3, sum); M5: 1-2 (reps 1,2); M6: 1",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"exclude" : ["Microglia","NPCs","OPCs","Astrocytes","Oligodendrocytes"],
	"outputPath": "runs/",
	"outputName": "Cortex_Neurons_Only.html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 5,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"replicate_mode": False,
	"excludeMarkers" : True,
	"minify": False,
}
results = t.TAP(**parameters)
