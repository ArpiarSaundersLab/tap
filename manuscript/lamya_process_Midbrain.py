import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from datetime import datetime

#open the data
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Midbrain.h5ad")
print(adata.shape)

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

#remove replicate n
adata = adata[adata.obs["replicates"] != "replicate_2", :]
adata.obs.replicates.value_counts()

parameters = {
	"adataObject": adata,
	"name": "Midbrain - All Replicated Combined",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicate","Cell Type"],
	"exclude" : ["Microglia","OPCs"], 
	"outputPath": "runs/",
	"outputName": "Lamya_Midbrain_Combined_"+str(datetime.now().strftime("%m-%d-%Y"))+".html",
	"excludeMarkers" : True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler", #randomoversampler or smote
	"minify": False,
	"removeOutliers": False,
	"replicate_mode": True,
}

results = t.TAP(**parameters)