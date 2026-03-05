import tap as t
import scanpy as sc
from datetime import datetime

#open the data
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad")

#tak a look at the value counts
adata.obs.replicates.value_counts()
adata.obs.cell_type.value_counts()

#replicates and cell types
parameters = {
	"adataObject": adata,
	"name": "Cortex - All Replicated Combined",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	"exclude" : ["Microglia","NPCs","OPCs"],
	"outputPath": "runs/",
	"outputName": "Lamya_Cortex_Combined_"+str(datetime.now().strftime("%m-%d-%Y"))+".html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler", #randomoversampler or smote
	"minify": False,
	"replicate_mode": True,
	"excludeMarkers" : True,
}
results = t.TAP(**parameters)