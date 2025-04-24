import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

#open the h5ad file
adata = sc.read_h5ad("../archive/data/GSE249416/GSE249416_AAV_ctxobj.Robj.h5ad")
adata.obs["replicate"] = np.random.randint(0, 2, size=adata.obs.shape[0])
adata.obs["replicate"] = adata.obs["replicate"].astype(str)
adata.obs["replicate"] = "replicate_" + adata.obs["replicate"]
adata.obs.groupby("replicate")["xincelltype230416"].value_counts()

parameters = {
	#"filename" : "../archive/data/GSE249416/GSE249416_AAV_ctxobj.Robj.h5ad",
	"adataObject": adata,
	"name": "Xin Jin 2024 - Test Replicate Mode",
	# "genes": ["BC1a","BC1b","BC1c",
	# 		  "BC2a","BC2b","BC2c",
	# 		  "BC3a","BC3b","BC3c",
	# 		  "BC4a","BC4b","BC4c",
	# 		  "BC5a","BC5b","BC5c",
	# 		  "BC6a","BC6b","BC6c",
	# 		  "BC7a","BC7b","BC7c",
	# 		  "BC8a","BC8b","BC8c",
	# 		  "BC9a","BC9b","BC9c",
	# 		  "BC10a","BC10b","BC10c",
	# 		  "BC11a","BC11b","BC11c",
	# 		  "BC12a","BC12b","BC12c",
	# 		  "BC13a","BC13b","BC13c",
	# 		  "BC14a","BC14b","BC14c",],
	"genes": ["BC1a","BC1b",],
	"exclude" : ["Microglia","Mural","Fibroblast"],
	"categories": ["replicate","xincelltype230416"],
	"categoryNames": ["Replicate","Cell Type"],
	"outputPath": ".",
	"outputName": "xin_test2.html",
	"excludeMarkers" : True,
	"mapOnly" : False,
	"showRF": True,
	"removeOutliers": True,
	"clusterMethod": "threshold",
	"clusterThreshold": 1,
	"clusterThresholdGreaterThanOrEqual": 2,
	"clusterThresholdLessThanOrEqual": 0,
	"rfType": "classifier",
	"rfHyperParameterTune": False,
	"rfHyperParameterIterations": 20,
	"rfPermuteFeatureImportance": False,
	"rfPermuteRepeats": 5,
	"balance": "smote",
	"minify": False,
}

results = t.TAP(**parameters)




#results.data_structure["replicate_0"]["ULPN"]
#results.data_structure["BC1a"]["replicate_1"]["ULPN"]


#Check if the RF is empty
def safe_gt(val, threshold):
    try:
        return float(val) > threshold
    except (ValueError, TypeError):
        return False


#iterate all nodes of the tree
for node in results.data_structure:

	print("========================🌶️")
	print(node)
	print("========================🌶️")

	#COLUMNS
	#are there any other nodes at this level with the same name?
	duplicates_columns = [n for n in results.data_structure if n == node]
	if len(duplicates_columns) > 1:
		print(f"\tDuplicates: {duplicates_columns}")
	

	#ROWS
	#are there any nodes at this level with the same name minus the last character?
	duplicates_rows = [n for n in results.data_structure if n[:-1] == node[:-1]]
	if len(duplicates_rows) > 1:

		#check the node.RF of the duplicates to see if any of the keys are the same and values  > 0 
		duplicates_rf = [list(filter(lambda x: results.data_structure[n]["RF"][x] > 0, results.data_structure[n]["RF"].keys())) for n in duplicates_rows]		
		duplicates_rf = set.intersection(*map(set, duplicates_rf))
		
		#check the node.DGE of the duplicate, but just look for the keys. They don't need to be > 0
		duplicates_dge = [list(results.data_structure[n]["DGE"].keys()) for n in duplicates_rows]
		duplicates_dge = set.intersection(*map(set, duplicates_dge))
		
		#which appear in both RF and DGE
		duplicates_both = [x for x in duplicates_rf if x in duplicates_dge]

		print(f"\tDuplicate rows: {duplicates_rows}")
		print(f"\tRF Row peppers:{duplicates_rf}")
		print(f"\tDGE Row peppers:{duplicates_dge}")
		print(f"\tBoth Row peppers:{duplicates_both}\n")
	

	#show all the keys this node that are in results.data_structure.keys()
	keys = [k for k in results.data_structure.keys() if k in results.data_structure[node].keys()]
	
	#iterate all of the node keys
	for k in keys:
		list_buffer_rf = []
		list_buffer_dge = []
		for row in duplicates_rows:
			list_buffer_rf.append([
				x for x in results.data_structure[row][k]['RF']
				if safe_gt(results.data_structure[row][k]['RF'][x], 0)
			])

			#same for DGE, but doesn't need to be there
			list_buffer_dge.append([
				x for x in results.data_structure[row][k]['DGE']
			])
			
		
		#intersection of all values in each list
		list_buffer_rf = list(set.intersection(*map(set, list_buffer_rf)))
		list_buffer_dge = list(set.intersection(*map(set, list_buffer_dge)))
		
		#which appear in both RF and DGE
		list_buffer_both = [x for x in list_buffer_rf if x in list_buffer_dge]

		print(f"\tRF {k}:{duplicates_rows}: peppers: {list_buffer_rf}")
		print(f"\tDGE {k}:{duplicates_rows}: peppers: {list_buffer_dge}")
		print(f"\tBoth {k}:{duplicates_rows}: peppers: {list_buffer_both}\n")




#1 - appears in RF across all replicates 
#2 - appears in RF and DGE across all replicates
#3 - The features are known receptors, exist in the plasma membrane, or extracellular. 