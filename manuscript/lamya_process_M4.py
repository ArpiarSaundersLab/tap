import tap as t
import scanpy as sc
from datetime import datetime

#open the data
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/molotAAV/Code/Analyses/M4/Cortex_1000_Kash_MT11/Data/Dialout/preprocessed.h5ad")

#take a look at the value counts
adata.obs.cell_type.value_counts()

#Cortex M6
parameters = {
	"adataObject": adata,
	"name": "Cortex: M4",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Retro_Kash","Retro_10","Rh10","Anc80","SHH10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"exclude" : ["Microglia","NPCs","OPCs"],
	"outputPath": "runs/",
	"outputName": "Lamya_Cortex_M4_"+str(datetime.now().strftime("%d-%m-%Y"))+".html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler", #randomoversampler or smote
	"replicate_mode": False,
	"excludeMarkers" : True,
}
results = t.TAP(**parameters)