import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tap as t
from memory_profiler import memory_usage
import os
import gc
import sys
import time
sns.set(style="whitegrid")
sc.settings.set_figure_params(dpi=300)

#opens the Cortex dataset
filename = "/home/kenny/Documents/OHSU/Projects/TAP/manuscript/data/Cortex.h5ad"

#defines the sample sizes to test
sizes = [2000, 5000, 10000, 20000]

#function for memory usage
def memory_test_tap(adata, *args, **kwargs):
	#params for the TAP run
	parameters = {
		"adataObject": adata,
		"name": "Runtime Tests",
		"genes": ["AAV1","AAV2","AAV5","AAV6","AAV7"],
		"categories": ["cell_type"],
		"categoryNames": ["Cell Type"],
		"exclude" : ["Microglia","NPCs","OPCs","Astrocytes","Oligodendrocytes"],
		"outputPath": "runs/",
		"outputName": "Runtime.html",
		"clusterMethod": "threshold",
		"clusterThreshold": 1,
		"balance": "randomoversampler",
		"replicate_mode": False,
		"excludeMarkers" : True,
	}

	#builds TAP
	results = t.TAP(**parameters)

	#clears the memory
	del results
	gc.collect()


#creates csv for results
df = pd.DataFrame(columns=["Size", "Iteration", "Runtime_Seconds","Memory_Usage"])

#iterates the sizes
for size in sizes:

	#read in the fresh data
	adata = sc.read_h5ad(filename)

	#determine highly variable genes
	sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)

	#downsample the dataset to 10k cells to reduce runtime
	sc.pp.sample(adata, n=size)

	#start time of the loop below
	for i in range(5):

		#set the start time 
		start_time = time.time()

		#run the memory test
		mem_usage = memory_usage((memory_test_tap, (adata,)))
		total_mem_usage = max(mem_usage) - min(mem_usage)

		#calcs the runtime
		runtime = time.time() - start_time
		
		#concat to dataframe
		df = pd.concat([df, pd.DataFrame({"Size":[size], 
										"Iteration":[i+1], 
										"Runtime_Seconds":[runtime]
										,"Memory_Usage":[total_mem_usage]})], 
										ignore_index=True)
		
		#writes to csv after each run
		df.to_csv("data/runtime_mem_results.csv", index=False)

		print(f"Size: {size} | Iteration: {i+1} | Runtime: {runtime:.2f} seconds | Memory Usage: {total_mem_usage:.2f} MiB")

		#clears memory
		gc.collect()



#opens the results csv
df = pd.read_csv("data/runtime_mem_results.csv")

#plots runtime for each size
plt.figure(figsize=(8,8))
sns.barplot(x="Size", y="Runtime_Seconds", data=df, palette="binary")
plt.title("TAP Runtime by Dataset Size\n", fontsize=30)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.xlabel("\nDataset Size (Cells)", fontsize=30)
plt.ylabel("Runtime (Seconds)\n", fontsize=30)
plt.xticks(rotation=90)
plt.show()

#plots memory usage for each size
plt.figure(figsize=(8,8))
sns.barplot(x="Size", y="Memory_Usage", data=df, palette="binary")
plt.title("TAP Memory Usage by Dataset Size\n", fontsize=30)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.xlabel("\nDataset Size (Cells)", fontsize=30)
plt.ylabel("Memory Usage(Mb)\n", fontsize=30)
plt.xticks(rotation=90)
plt.show()
