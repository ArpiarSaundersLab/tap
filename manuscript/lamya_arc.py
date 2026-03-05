import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from datetime import datetime

#######################################################################################
# Run Cortex - Markers -False
#######################################################################################
#open the data
adata = sc.read_h5ad("data/Cortex.h5ad")

#replicates and cell types
parameters = {
	"adataObject": adata,
	"name": "Cortex - All Replicates",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	"exclude" : ["Microglia"],
	"outputPath": "runs/",
	"outputName": "Lamya_Cortex_IncludeMarkers_"+str(datetime.now().strftime("%m-%d-%Y"))+".html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler", #randomoversampler or smote
	"minify": False,
	"replicate_mode": True,
	"excludeMarkers" : False,
}
results = t.TAP(**parameters)


#######################################################################################
# Run Cortex  - Markers -True
#######################################################################################
#open the data
adata = sc.read_h5ad("data/Cortex.h5ad")

#replicates and cell types
parameters = {
	"adataObject": adata,
	"name": "Cortex - All Replicates",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Retro","Rh10","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	"exclude" : ["Microglia"],
	"outputPath": "runs/",
	"outputName": "Lamya_Cortex_ExcludeMarkers_"+str(datetime.now().strftime("%m-%d-%Y"))+".html",
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



#######################################################################################
# Run Midbrain  - Markers -True
#######################################################################################

#open the data
adata = sc.read_h5ad("data/Midbrain.h5ad")

#tak a look at the value counts
adata.obs.replicates.value_counts()
adata.obs.cell_type.value_counts()

#rename OPCs/Oligos to Oligos and OPCs to Oligos (markers are indicate either)
adata.obs["cell_type"] = adata.obs["cell_type"].replace({"OPCs/Oligos":"Oligos"})
adata.obs["cell_type"] = adata.obs["cell_type"].replace({"OPCs":"Oligos"})

#rename NPC/Excitatory to Excitatory(markers indicate either)
adata.obs["cell_type"] = adata.obs["cell_type"].replace({"NPC/Excitatory":"Excitatory"})

#take a look at cell type counts again
adata.obs.cell_type.value_counts()

#remove replicate 5
adata = adata[adata.obs["replicates"] != "replicate_5", :]

parameters = {
	"adataObject": adata,
	"name": "Midbrain - All Replicates",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	"exclude" : ["Microglia"], 
	"outputPath": "runs/",
	"outputName": "Lamya_Midbrain_ExcludeMarkers"+str(datetime.now().strftime("%m-%d-%Y"))+".html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler", #randomoversampler or smote
	"minify": False,
	"removeOutliers": False,
	"replicate_mode": True,
	"excludeMarkers" : True,
}

results = t.TAP(**parameters)

#######################################################################################
# Run Midbrain  - Markers -False
#######################################################################################



#open the data
adata = sc.read_h5ad("data/Midbrain.h5ad")

#tak a look at the value counts
adata.obs.replicates.value_counts()
adata.obs.cell_type.value_counts()

#rename OPCs/Oligos to Oligos and OPCs to Oligos (markers are indicate either)
adata.obs["cell_type"] = adata.obs["cell_type"].replace({"OPCs/Oligos":"Oligos"})
adata.obs["cell_type"] = adata.obs["cell_type"].replace({"OPCs":"Oligos"})

#rename NPC/Excitatory to Excitatory(markers indicate either)
adata.obs["cell_type"] = adata.obs["cell_type"].replace({"NPC/Excitatory":"Excitatory"})

#take a look at cell type counts again
adata.obs.cell_type.value_counts()

#remove replicate 5
adata = adata[adata.obs["replicates"] != "replicate_5", :]

parameters = {
	"adataObject": adata,
	"name": "Midbrain - All Replicates",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"excludeGenes": ["Retro_Kash","Retro_10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	"exclude" : ["Microglia"], 
	"outputPath": "runs/",
	"outputName": "Lamya_Midbrain_IncludeMarkers"+str(datetime.now().strftime("%m-%d-%Y"))+".html",
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler", #randomoversampler or smote
	"minify": False,
	"removeOutliers": False,
	"replicate_mode": True,
	"excludeMarkers" : False,
}

results = t.TAP(**parameters)
