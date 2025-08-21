import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad


adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Midbrain.h5ad")
adata = adata[~adata.obs["replicates"].isin(["replicate_2"]),:]

parameters = {
	"adataObject": adata,
	"name": "Midbrain: M4: 1; M5: 2-3; M6: 1 = reps 1-5",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"categories": ["replicates","cell_type"],
	"categoryNames": ["Replicates","Cell Type"],
	"exclude" : ["Microglia","NPCs","OPCs","NPC/Excitatory"], #m5rep1 no signal
	"outputPath": "runs/",
	"outputName": "Midbrain_M1-6-2.html",
	"excludeMarkers" : True,
	"clusterMethod": "threshold", 
	"clusterThreshold": 3,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "randomoversampler",
	"minify": False,
	"replicate_mode": True,
	"removeOutliers": False,
}

results = t.TAP(**parameters)


adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Midbrain.h5ad")

parameters = {
	"adataObject": adata,
	"name": "Midbrain - All Replicated Combined",
	"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7","AAV8","AAV9","Rh10","Retro","Anc80","SHH10"],
	"categories": ["cell_type"],
	"categoryNames": ["Cell Type"],
	"exclude" : ["Microglia","NPCs","OPCs"], 
	"outputPath": "runs/",
	"outputName": "Midbrain_Combined.html",
	"excludeMarkers" : True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 3,
	"clusterThresholdLessThanOrEqual": 0,
	"balance": "smote",
	"minify": False,
	"replicate_mode": True,
	"removeOutliers": False,
}

results = t.TAP(**parameters)