import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import random
import json
import tables
import copy
import gc
import os
import anndata as ad
import warnings
import time
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.sparse import issparse
from scipy import stats
import scipy.sparse as sp
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.feature_selection import SelectFromModel
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GridSearchCV
from sklearn.utils import shuffle
from sklearn.mixture import GaussianMixture
from matplotlib.colors import ListedColormap
import pkg_resources
import shutil
from tap.FeatureSelection import FeatureSelection
warnings.filterwarnings("ignore")
sns.set(style="whitegrid")
sc.settings.verbosity = 0
from tqdm import tqdm
from biothings_client import get_client

class TAP:

	colors = [(0.8, 0.8, 0.8)] + [(plt.cm.Reds(i / 255)) for i in range(80, 256)]
	cmap_custom = ListedColormap(colors)

	def __init__(self, name="Heatmap", filename=None, adataObject=None, categories=[], categoryNames=[], 
	genes=[], exclude=[], useRaw=False, outputPath=".", outputName="tap.html", cpus = 2, 
	mapOnly=False, showDE=True, showRF=True, showUMAP=True, showDetails=True, 
	clusterMethod="threshold", clusterThreshold=1, useLog10=False, showPlots=False, 
	excludeGenes=[], deMethod="wilcoxon", rfHyperParameterTune=False, rfHyperParameterIterations=50, rfType="classifier",
	balance=None, useAllGenes=False, rfPermuteFeatureImportance=False, rfPermuteRepeats=2,
	runCellTypist=False, cellTypistModel="Mouse_Whole_Brain.pkl", cellTypistPlots=False, cellTypistLevels=3,
	minify = True, excludeMarkers=False, removeOutliers=True, minCells=20, minSeroTypeCells=20,
	clusterThresholdGreaterThanOrEqual=None, clusterThresholdLessThanOrEqual=None):
		self.name = name 
		self.filename = filename
		self.adataObject = adataObject
		self.categories = categories 
		self.exclude = exclude 
		self.genes = genes
		self.result = None
		self.df = None
		self.accuracy = None
		self.outputPath = outputPath
		self.outputName = outputName
		self.balance = balance
		self.useAllGenes = useAllGenes
		self.rfPermuteFeatureImportance = rfPermuteFeatureImportance
		self.rfPermuteRepeats = rfPermuteRepeats
		self.rfType = rfType
		self.removeOutliers = removeOutliers
		self.minCells = minCells
		self.minSeroTypeCells = minSeroTypeCells
		self.totalIterations = 0
		if self.outputPath[-1] != "/":
			self.outputPath = self.outputPath + "/"
		if not os.path.exists(self.outputPath):
			os.makedirs(self.outputPath)
		self.cpus = cpus
		self.mapOnly = mapOnly
		self.categoryNames = categoryNames
		self.showDE = showDE
		self.showRF = showRF
		self.showUMAP = showUMAP
		self.showDetails = showDetails
		self.clusterMethod = clusterMethod
		self.clusterThreshold = clusterThreshold
		self.clusterThresholdGreaterThanOrEqual = clusterThresholdGreaterThanOrEqual
		self.clusterThresholdLessThanOrEqual = clusterThresholdLessThanOrEqual
		self.excludeGenes = excludeGenes
		self.useLog10 = useLog10
		self.deMethod = deMethod
		self.showPlots = showPlots
		self.rfHyperParameterTune = rfHyperParameterTune
		self.rfHyperParameterIterations = rfHyperParameterIterations
		self.runCellTypist = runCellTypist
		self.cellTypistModel = cellTypistModel
		self.cellTypistPlots = cellTypistPlots
		self.cellTypistLevels = cellTypistLevels
		self.excludeMarkers = excludeMarkers
		self.markers ={}
		self.copyTapTemplate()
		self.loadData()
		
		if self.runCellTypist == True:
			self.generateCellTypistPredictions()
		
		self.runTimeEstimate()

		self.generateHeatmapMetadata()

		if self.excludeMarkers == True:
			self.determineMarkers()
		if self.mapOnly == True:
			self.generateCSV()
		else:
			self.generateCSV()
			self.generateJSON()
		if minify == True:	
			self.minifyHtml()

	def copyTapTemplate(self):		
		shutil.copy(pkg_resources.resource_filename("tap", "templates")+"/tap.html", self.outputPath + self.outputName)



	def loadData(self):

		#if user passes in a specific object name, use that object
		if self.adataObject:
			self.adata = self.adataObject
		else:
			self.adata = sc.read(self.filename)

		#does the raw object exist?
		if self.adata.raw == None:
			raise ValueError("Error: No raw.X found in the AnnData object.")


		#is this a sparse matrix?
		if issparse(self.adata.X) == True:
			self.adata.X = self.adata.X.toarray()
			#self.adata.raw = ad.AnnData(self.adata.raw.X.toarray())
			self.adata.raw = ad.AnnData(
				X=self.adata.raw.X.toarray(),
				var=self.adata.raw.var.copy(),
			)


		#check if the raw object has float values. If so, they aren't raw counts.
		if self.adata.raw.X[0].dtype == "float32" or self.adata.raw.X[0].dtype == "float64":
			if np.all(self.adata.raw.X[0] % 1 != 0):
				raise ValueError("Error: raw.X appears to be normalized or log transformed. Please use raw.X.")

		input_adata_genes = self.adata.var_names.tolist()

		#create a copy of the adata object
		self.adata_bak = self.adata.copy()
		
		#convert to raw
		self.adata = self.adata.raw.to_adata() 

		#we only want the original subset of genes, but the raw values of those genes
		if self.useAllGenes == False:
			
			#check if the input_adata_genes are in the adata
			input_adata_genes = [x for x in input_adata_genes if x in self.adata.var_names.tolist()]
			
			#select the genes
			self.adata = self.adata[:, input_adata_genes]

		#iterate each gene and check if it exists in the adata object, then add from obs if necessary
		for gene in self.genes:
			if gene not in self.adata.var_names.tolist() and gene in self.adata.obs.columns.tolist():
				print(f"Notice: {gene} not found in adata object, but in obs. The counts should be raw counts. Adding...")
				df = pd.DataFrame({"Cell": self.adata.obs.index.tolist(), f"{gene}": self.adata.obs[gene].values})
				df = df.set_index("Cell")
				adata2 = ad.AnnData(df)
				self.adata = ad.concat([self.adata, adata2], join="outer", axis=1, merge="only")
				del adata2
				gc.collect()

		#make var names unique
		self.adata.var_names_make_unique()

		#calculate the percent of counts from mitochondrial genes
		self.adata.var["mt"] = self.adata.var_names.str.startswith("mt-")
		sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
		#let users decide ahead of time
		#self.adata = self.adata[self.adata.obs.pct_counts_mt < 10, :]

		#calculate the infection_count
		self.adata.obs['infection_count'] = np.sum(self.adata[:, self.genes].to_df(), axis=1).astype(int)
		is_infected = self.adata.obs['infection_count'] > 0
		self.adata.obs['infected'] = pd.Categorical(is_infected.astype(str))

	def runTimeEstimate(self):
		category_length = len(self.categories)
		for category in self.adata.obs[self.categories[0]].unique():
			if category not in self.exclude:
				self.totalIterations +=1
				if category_length > 1:
					for subcategory in self.adata[self.adata.obs[self.categories[0]] == category, :].obs[self.categories[1]].unique():
						self.totalIterations +=1
						if category_length > 2:
							for subsubcategory in self.adata[(self.adata.obs[self.categories[0]] == category) & (self.adata.obs[self.categories[1]] == subcategory), :].obs[self.categories[2]].unique():
								self.totalIterations +=1
		self.totalIterations = self.totalIterations+(self.totalIterations*len(self.genes)) + len(self.genes) + 1

	def generateCellTypistPredictions(self, showPlots=False):
		
		#import the celltypist libraries
		import celltypist
		from celltypist import models

		#set the model path
		models.models_path = pkg_resources.resource_filename("tap", "models")

		#download relevant model
		models.download_models(model = self.cellTypistModel)

		#deep copy the adata object
		self.celltypist_adata = copy.deepcopy(self.adata_bak)

		#load the raw object to a celltypist object
		self.celltypist_adata = self.celltypist_adata.raw.to_adata() 

		#CT expects normalized to 10,000 and natural log scale
		sc.pp.normalize_total(self.celltypist_adata, target_sum=1e4)
		sc.pp.log1p(self.celltypist_adata)

		#does uns neighbors exist?
		if "neighbors" not in self.celltypist_adata.uns.keys():
			sc.pp.neighbors(self.celltypist_adata)

		#run celltypist on the data
		predictions = celltypist.annotate(self.celltypist_adata, model = self.cellTypistModel,
										majority_voting = True, 
										mode="best match")
		
		#save the predictions
		self.celltypist_adata = predictions.to_adata()

		#plot if desired
		if self.cellTypistPlots == True:
			
			#plot the dotplot
			celltypist.dotplot(predictions, 
								use_as_reference =self.categories[0], 
								use_as_prediction = 'majority_voting', 
								title='CellTypist Predictions')
			
			#take a look at the confidence scores
			self.celltypist_adata.obs.conf_score.hist()

			#plot umap
			sc.pl.umap(self.celltypist_adata, color = [self.categories[0], 'majority_voting'])

		#add obs placeholders to the molotAAV object
		self.celltypist_adata.obs["celltypist_primary"] = ""
		self.celltypist_adata.obs["celltypist_secondary"] = ""
		self.celltypist_adata.obs["celltypist_tertiary"] = ""

		#iterate through the celltypist object and add the celltypist predictions
		if self.cellTypistModel == "Mouse_Whole_Brain.pkl":
			
			#open a csv as a pandas dataframe
			abc_atlas_metadata = pd.read_csv(pkg_resources.resource_filename("tap", "models")+"/abc_atlas_metadata.csv")

			#iterate each cell in the object and query the subclass_id_label column
			for i in range(len(self.celltypist_adata.obs)):
				majority_voting = self.celltypist_adata.obs["majority_voting"][i]
				if majority_voting in abc_atlas_metadata["subclass_id_label"].values:
					abc_subclass_df = abc_atlas_metadata[abc_atlas_metadata["subclass_id_label"] == majority_voting]
					self.celltypist_adata.obs["celltypist_primary"][i] = abc_subclass_df["nt_type_label"].values[0]
					self.celltypist_adata.obs["celltypist_secondary"][i] = abc_subclass_df["class_label"].values[0]
					self.celltypist_adata.obs["celltypist_tertiary"][i] = abc_subclass_df["subclass_label"].values[0]

			#save the celltypist object obs to the adata object
			self.adata.obs["celltypist_primary"] = self.celltypist_adata.obs["celltypist_primary"].str.replace(":","").str.replace("-","")
			self.adata.obs["celltypist_secondary"] = self.celltypist_adata.obs["celltypist_secondary"].str.replace(":","").str.replace("-","")
			self.adata.obs["celltypist_tertiary"] = self.celltypist_adata.obs["celltypist_tertiary"].str.replace(":","").str.replace("-","")
			
			if self.cellTypistLevels == 1:
				self.categories = ["celltypist_primary"]
				self.categoryNames = ["NT Type"]
			elif self.cellTypistLevels == 2:
				self.categories = ["celltypist_primary","celltypist_secondary"]
				self.categoryNames = ["NT Type","Class"]
			elif self.cellTypistLevels == 3:
				self.categories = ["celltypist_primary","celltypist_secondary","celltypist_tertiary"]
				self.categoryNames = ["NT Type","Class","Subclass"]

			#replace any None values from the celltypist values above with "NA"
			self.adata.obs["celltypist_primary"] = self.adata.obs["celltypist_primary"].replace("", "NA")
			self.adata.obs["celltypist_secondary"] =self.adata.obs["celltypist_secondary"].replace("", "NA")
			self.adata.obs["celltypist_tertiary"] = self.adata.obs["celltypist_tertiary"].replace("", "NA")
			self.adata.obs["celltypist_primary"] = pd.Categorical(self.adata.obs["celltypist_primary"].fillna("NA"))
			self.adata.obs["celltypist_secondary"] = pd.Categorical(self.adata.obs["celltypist_secondary"].fillna("NA"))
			self.adata.obs["celltypist_tertiary"] = pd.Categorical(self.adata.obs["celltypist_tertiary"].fillna("NA"))

		else:
			for i in range(len(self.celltypist_adata.obs)):
				majority_voting = self.celltypist_adata.obs["majority_voting"][i]
				self.celltypist_adata.obs["celltypist_primary"][i] = self.celltypist_adata.obs["majority_voting"][i]

			#save the celltypist object obs to the adata object
			self.adata.obs["celltypist_primary"] = self.celltypist_adata.obs["celltypist_primary"].str.replace(":","")
			self.adata.obs["celltypist_primary"] = self.adata.obs["celltypist_primary"].replace("", "NA")
			self.adata.obs["celltypist_primary"] = pd.Categorical(self.adata.obs["celltypist_primary"].fillna("NA"))
			self.categories = ["celltypist_primary"]
			self.categoryNames = ["CellTypist"]

		#delete the celltypist object
		del self.celltypist_adata
		gc.collect()


	def generateHeatmapMetadata(self):
		metadata = {
			"filename": self.filename,
			"title": self.name,
			"outputName": self.outputName,
			"categories": self.categories,
			"categoryNames": self.categoryNames if self.categoryNames else self.categories,			
			"genes": self.genes,
			"DE": self.showDE,
			"RF": self.showRF,
			"showPlots": self.showPlots,
			"UMAP": self.showUMAP,
			"Details": self.showDetails,
			"clusterMethod": self.clusterMethod,
			"clusterThreshold": self.clusterThreshold,
			"clusterThresholdGreaterThanOrEqual": self.clusterThresholdGreaterThanOrEqual,
			"clusterThresholdLessThanOrEqual": self.clusterThresholdLessThanOrEqual,
			"deMethod": self.deMethod,
			"mapOnly": self.mapOnly,
			"excludeMarkers": self.excludeMarkers,
			"cellTypist": self.runCellTypist,
			"cellTypistModel": self.cellTypistModel,
			"cellTypistLevels": self.cellTypistLevels,
			"rfType": self.rfType,
			"rfHyperParameterTune": self.rfHyperParameterTune,
			"rfHyperParameterIterations": self.rfHyperParameterIterations,
			"permute": self.rfPermuteFeatureImportance,
			"permuteRepeats": self.rfPermuteRepeats,
		}

		#instead of writing to file, open the tap html file and replace the METADATA placeholder with proper json
		with open(self.outputPath+self.outputName, 'r') as file:
			filedata = file.read()

		#replace the target string
		filedata = filedata.replace('METADATA', json.dumps(metadata))
		if self.mapOnly == True:
			filedata = filedata.replace('CLUSTERING_DETAILS', "''")

		#write the file out again
		with open(self.outputPath+self.outputName, 'w') as file:
			file.write(filedata)

		#line break for each key value pair
		with open(f'{self.outputPath}metadata_{self.clean_string(self.name)}.txt', 'w') as file:
			file.write(json.dumps(metadata, indent=4))

	def clean_string(self, string):
		import re
		return re.sub(r"[/\\?%*:|\"<>\x7F\x00-\x1F]", "-",string)


	def determineMarkers(self):
		adata_copy = copy.deepcopy(self.adata)

		for category in self.categories:
			try:
				sc.tl.rank_genes_groups(adata_copy, 
										groupby=category,
										key_added=category, 
										method=self.deMethod)
			except Exception as e:
				print(f"Error in determineMarkers: {e}")
		
		for category in self.categories:
			categories = adata_copy.obs[category].unique()
			try:
				for catName in categories:
					self.markers[catName] = sc.get.rank_genes_groups_df(adata_copy, 
																		key=category, 
																		group=catName)
					#keep only the markers with a score of 2.2 std deviations above the mean.
					self.markers[catName] = self.markers[catName][self.markers[catName]["scores"] > self.markers[catName]["scores"].mean() + 2.2*self.markers[catName]["scores"].std()]
					self.markers[catName] = self.markers[catName][self.markers[catName]["pvals_adj"] <= 0.05]
					self.markers[catName] = self.markers[catName].sort_values(by="scores", ascending=False)
					self.markers[catName] = self.markers[catName].reset_index(drop=True)
			except Exception as e:
				print(f"Error in determineMarkers2: {e}")
		
		del adata_copy
		gc.collect()


	def generateCSV(self, plot=False):

		# make a deep copy of self.result object
		self.result_copy = copy.deepcopy(self.adata)

		#normalize by library size using Scanpy's normalize_total function
		sc.pp.normalize_total(self.adata, target_sum=1e4)
		
		#calculate the infection_count again after normalization
		self.adata.obs['infection_count'] = np.sum(self.adata[:, self.genes].to_df(), axis=1).astype(int)
		is_infected = self.adata.obs['infection_count'] > 0
		self.adata.obs['infected'] = is_infected.astype(str)

		#remove any primary cats in the exclude list
		self.adata = self.adata[~self.adata.obs[self.categories[0]].isin(self.exclude), :]
		
		#if there's a secondary category or tertiary, remove those too
		if len(self.categories) > 1:
			self.adata = self.adata[~self.adata.obs[self.categories[1]].isin(self.exclude), :]
		
		if len(self.categories) > 2:
			self.adata = self.adata[~self.adata.obs[self.categories[2]].isin(self.exclude), :]

		#only keep the genes
		self.adata = self.adata[:,self.genes]
		
		#remove any outliers from the infection_count obs by looking at the quantiles
		if self.removeOutliers == True:
			self.adata = self.adata[self.adata.obs["total_counts"] < self.adata.obs["total_counts"].quantile(0.98), :]
			self.adata = self.adata[self.adata.obs["infection_count"] < self.adata.obs["infection_count"].quantile(0.98), :]

		#plot a histogram of the infection_count 
		if plot == True:
			self.adata.obs["total_counts"].hist(bins=100)
			self.adata.obs["infection_count"].hist(bins=100)

		#scale the infection_count data with a max of 1
		max_infection_count = self.adata.obs['infection_count'].max()
		self.adata.obs['infection_count_normalized'] = (
			self.adata.obs['infection_count'] / max_infection_count
		)

		#if sparse, convert to dense. 
		if issparse(self.adata.X) == True:
			self.adata.X = self.adata.X.toarray()
			self.adata.X = self.adata.X.astype(float)

		#iterate all and find the max value across all the genes
		max_value = np.max([
			self.adata[:, gene].X.max()
			for gene in self.genes
		]).astype(float)

		#normalize each of the serotypes the same way
		for serotype in self.genes:

			#extract the column for the serotype
			serotype_data = self.adata[:, serotype].X.toarray().flatten()

			#find the maximum value in the serotype data
			#max_value = serotype_data.max().astype(float)
			
			#normalize if max_value is not zero... avoid division by zero
			if max_value != 0:

				#normalize the gene data
				normalized_data = (serotype_data / max_value)

				#reshape normalized_data to a 2D array with one column
				normalized_data = normalized_data.reshape(-1, 1)

				#assign the normalized_data to the serotype column
				self.adata[:, serotype].X = normalized_data

		#create a df with the neccessary columns
		df = pd.DataFrame(self.adata.obs[self.categories+['infection_count','total_counts','infection_count_normalized']])
		df.index.name = "Cell"

		#add the serotype data to the df
		df = pd.concat([df, self.adata.to_df()[self.genes]], axis=1)

		#get the cell umap coordinates
		df["umap_x"] = self.adata.obsm['X_umap'][:,0]
		df["umap_y"] = self.adata.obsm['X_umap'][:,1]

		#create a list of True values for each category
		category_bools = []
		for category in self.categories:
			category_bools.append(True)

		#order by highest to lowest; then, infection_count.
		df = df.sort_values(by=self.categories+['infection_count'],
							ascending=category_bools+[False])
		
		#save the df
		self.df = df

		#rather than saving the csv, write the csv to the tap.html file with line breaks as \n
		with open(self.outputPath+self.outputName, 'r') as file:
			filedata = file.read()

		#replace the target string
		filedata = filedata.replace('HEATMAP_DATA', df.to_csv(index=False,na_rep='NA').replace('\n', '\\n'))

		#write the file out again
		with open(self.outputPath+self.outputName, 'w') as file:
			file.write(filedata)

	def minifyHtml(self):
		import minify_html

		#open the tap.html file
		with open(self.outputPath+self.outputName, 'r') as file:
			filedata = file.read()
		
		#minify the html
		minified = minify_html.minify(filedata,
									minify_js=True,
									minify_css=True,
									keep_closing_tags=False,
									preserve_brace_template_syntax=False,
									keep_input_type_text_attr = False,
									)

		#replace only the text "\\n" with "\n"
		minified = minified.replace('\\\\n', '\\n')

		# Write the file out again
		with open(self.outputPath+self.outputName, 'w') as file:
			file.write(minified)


	def generateJSON(self):
		#run the GMM+RF/DE analysis on all the data
		adata_feature_selection = FeatureSelection(self.result_copy, 
								serotype_of_interest="infection_count", 
								outputPath=self.outputPath,
								primary=None, 
								min_cells=self.minCells, 
								min_serotype_cells=self.minSeroTypeCells, 
								threads=self.cpus,
								cluster_method=self.clusterMethod,
								cluster_threshold=self.clusterThreshold,
								clusterThresholdGreaterThanOrEqual=self.clusterThresholdGreaterThanOrEqual,
								clusterThresholdLessThanOrEqual=self.clusterThresholdLessThanOrEqual,
								showDE=self.showDE,
								showRF=self.showRF,
								showPlots=self.showPlots,
								category_observations=self.categories,
								exclude_category_names=self.exclude,
								serotype_list=self.genes,
								exclude_genes=self.excludeGenes,
								deMethod=self.deMethod,
								rfHyperParameterTune=self.rfHyperParameterTune,
								rfHyperParameterIterations = self.rfHyperParameterIterations,
								rfType=self.rfType,
								balance=self.balance,
								permute_feature_importance=self.rfPermuteFeatureImportance,
								permutation_repeats=self.rfPermuteRepeats,
								excludeMarkers=self.excludeMarkers,
								markers=self.markers,
								removeOutliers=self.removeOutliers,
								totalIterations=self.totalIterations,
								)

		#define a json structure
		data_structure = {}

		#define a structure ffor the base64 encoded images
		images = {}

		try:
			#dge data dictionary
			dge_all = adata_feature_selection.de_df.set_index('1_name')['0_name'].to_dict()
			dge_pval_adj = adata_feature_selection.de_df[['1_pval_adj','0_pval_adj']].values.tolist()
			dge_logfc = adata_feature_selection.de_df[['1_logfc','0_logfc']].values.tolist()

			#rf dictionary
			rf_all = adata_feature_selection.fi_df[:50].set_index('Feature')['Importance'].to_dict()

			#add the "All" key to the data_structure
			data_structure["All"] = {
				"RF": rf_all,
				"DGE": dge_all,
				"DGE_pval_adj": dge_pval_adj,
				"DGE_logfc": dge_logfc,
				"RF_A": adata_feature_selection.accuracy,
				"RF_AUC": adata_feature_selection.auc_score,
				"RF_BP": adata_feature_selection.rf_params,
			}
		except Exception as e:
			warnings.warn(f"Warning (0): {e}")
			dge_all = {"Error": "All"}
			rf_all = {"Error": "All"}
			dge_pval_adj = {"Error": "All"}
			dge_logfc = {"Error": "All"}
			adata_feature_selection.clustering_details = ""

		#image key for the "All" key
		images["All"] = adata_feature_selection.clustering_details

		#deep copy the adata_feature_selection object
		rf_obj = copy.deepcopy(adata_feature_selection)

		#run for each main primary, secondary, and tertiary category
		data_structure = self.generateRfDgeStructure(data_structure, node=False, rf_obj=rf_obj, images=images)

		#run for each gene (node)
		for gene in self.genes:
			currentIterationTemp = rf_obj.currentIteration
			rf_obj = copy.deepcopy(adata_feature_selection)
			rf_obj.currentIteration = currentIterationTemp
			data_structure = self.generateRfDgeStructure(data_structure, node=gene, rf_obj=rf_obj, images=images)
			gc.collect()

		#save to property
		self.data_structure = data_structure

		#rather than saving the json, write the json to the tap.html file as json on one line
		with open(self.outputPath+self.outputName, 'r') as file:
			filedata = file.read()
		
		#Replace the target string
		filedata = filedata.replace('FEATURE_DATA', json.dumps(data_structure))
		filedata = filedata.replace('CLUSTERING_DETAILS', json.dumps(images))
		
		#replace the word Infinity with 0 to play nice with javascript
		filedata = filedata.replace('Infinity', '0')


		#find the "peppers"
		print("Finding Peppers...")
		self.pepper_results = self.find_pepper_results()
		self.pepper_unique_list = self.get_unique_peppers(self.pepper_results)
		self.peppers = self.get_pepper_gene_details(self.pepper_unique_list)
		self.pepper_results2 = self.increase_pepper_heat(self.pepper_results, self.peppers)
		filedata = filedata.replace('PEPPER_DATA', json.dumps(self.pepper_results2))


		# Write the file out again
		with open(self.outputPath+self.outputName, 'w') as file:
			file.write(filedata)
		

		#line break for each key value pair
		with open(f'{self.outputPath}peppers_{self.clean_string(self.name)}.txt', 'w') as file:
			file.write(json.dumps(self.pepper_results2, indent=4))

		
		print(f"Finished Running")


	def generateRfDgeStructure(self, data_structure, node=False, rf_obj=False, images = {}):

		if self.useLog10 == True:
			gene_of_interest = "infection_count_log10"
		else:
			gene_of_interest = "infection_count"

		#add the node to the data_structure if it exists
		if node:

			try:
				
				node_lookup = node
				if self.useLog10 == True and node.lower() != "infection_count" and node.lower() != "infection_count_log10":
					node_lookup = node.lower() + "_log10"
				
				# change the focus of the rf_obj
				rf_obj.change_focus(primary=None,serotype_of_interest=node_lookup,min_cells=self.minCells,min_serotype_cells=self.minSeroTypeCells)

				# dge data dictionary
				dge_1 = rf_obj.de_df.set_index('1_name')['0_name'].to_dict()
				dge_1_pval_adj = rf_obj.de_df[['1_pval_adj','0_pval_adj']].values.tolist()
				dge_1_logfc = rf_obj.de_df[['1_logfc','0_logfc']].values.tolist()

				# rf dictionary
				rf_1 = rf_obj.fi_df[:50].set_index('Feature')['Importance'].to_dict()

			except Exception as e:
				warnings.warn(f"Warning(node level) (1): {node} {e}")
				dge_1 = {"Error": node}
				rf_1 = {"Error": node}
				dge_1_pval_adj = {"Error": node}
				dge_1_logfc = {"Error": node}
				rf_obj.clustering_details = ""


			# create the primary category key
			data_structure[node] = {
				"RF": rf_1,
				"DGE": dge_1,
				"DGE_pval_adj": dge_1_pval_adj,
				"DGE_logfc": dge_1_logfc,
				"RF_A": rf_obj.accuracy,
				"RF_AUC": rf_obj.auc_score,
				"RF_BP": rf_obj.rf_params,
			}

			images[str(node)] = rf_obj.clustering_details
			gene_of_interest = node


		#convert to lowercase
		gene_of_interest = gene_of_interest.lower()

		if self.useLog10 == True and gene_of_interest != "infection_count_log10":
			gene_of_interest = gene_of_interest + "_log10"

		#get the unique primaries
		primaries = self.df[self.categories[0]].unique()

		#iterate through the primaries
		for primary in primaries:

			try:

				# change the focus of the rf_obj
				rf_obj.change_focus(primary=primary,serotype_of_interest=gene_of_interest,min_cells=self.minCells,min_serotype_cells=self.minSeroTypeCells)

				# dge data dictionary
				dge_1 = rf_obj.de_df.set_index('1_name')['0_name'].to_dict()
				dge_1_pval_adj = rf_obj.de_df[['1_pval_adj','0_pval_adj']].values.tolist()
				dge_1_logfc = rf_obj.de_df[['1_logfc','0_logfc']].values.tolist()

				# rf dictionary
				rf_1 = rf_obj.fi_df[:50].set_index('Feature')['Importance'].to_dict()

			except Exception as e:
				warnings.warn(f"Warning (2): {node}:{primary}:{gene_of_interest} {e}")
				dge_1 = {"Error": primary}
				rf_1 = {"Error": primary}
				dge_1_pval_adj = {"Error": primary}
				dge_1_logfc = {"Error": primary}
				rf_obj.clustering_details = ""


			# create the primary category key
			if node:

				data_structure[node][primary] = {
					"RF": rf_1,
					"DGE": dge_1,
					"DGE_pval_adj": dge_1_pval_adj,
					"DGE_logfc": dge_1_logfc,
					"RF_A": rf_obj.accuracy,
					"RF_AUC": rf_obj.auc_score,
					"RF_BP": rf_obj.rf_params,
				}

				images[str(node)+'-'+str(primary)] = rf_obj.clustering_details

			else:

				data_structure[primary] = {
					"RF": rf_1,
					"DGE": dge_1,
					"DGE_pval_adj": dge_1_pval_adj,
					"DGE_logfc": dge_1_logfc,
					"RF_A": rf_obj.accuracy,
					"RF_AUC": rf_obj.auc_score,
					"RF_BP": rf_obj.rf_params,
				}

				images[str(primary)] = rf_obj.clustering_details

			# only run if there are more than 2 categories
			if len(self.categories) > 1:

				#get the unique secondaries for the cell type
				secondaries = self.df[self.df[self.categories[0]] == primary][self.categories[1]].unique()

				#iterate through the secondarys
				for secondary in secondaries:
				
					try:

						# change the focus of the rf_obj
						rf_obj.change_focus(primary=primary,secondary=secondary,serotype_of_interest=gene_of_interest,min_cells=self.minCells,min_serotype_cells=self.minSeroTypeCells)

						# dge data dictionary
						dge_2 = rf_obj.de_df.set_index('1_name')['0_name'].to_dict()
						dge_2_pval_adj = rf_obj.de_df[['1_pval_adj','0_pval_adj']].values.tolist()
						dge_2_logfc = rf_obj.de_df[['1_logfc','0_logfc']].values.tolist()


						# rf dictionary
						rf_2 = rf_obj.fi_df[:50].set_index('Feature')['Importance'].to_dict()

					except Exception as e:
						warnings.warn(f"Warning (3): {node} {primary}:{secondary}:{gene_of_interest} {e} ")
						dge_2 = {"Error": primary+":"+secondary}
						rf_2 = {"Error": primary+":"+secondary}
						dge_2_pval_adj = {"Error": primary+":"+secondary}
						dge_2_logfc = {"Error": primary+":"+secondary}
						rf_obj.clustering_details = ""


					# create the seconday category key
					if node:

						data_structure[node][primary][secondary] = {
							"RF": rf_2,
							"DGE": dge_2,
							"DGE_pval_adj": dge_2_pval_adj,
							"DGE_logfc": dge_2_logfc,
							"RF_A": rf_obj.accuracy,
							"RF_AUC": rf_obj.auc_score,
							"RF_BP": rf_obj.rf_params,
						}
						
						images[str(node)+'-'+str(primary)+'-'+str(secondary)] = rf_obj.clustering_details

					else:

						data_structure[primary][secondary] = {
							"RF": rf_2,
							"DGE": dge_2,
							"DGE_pval_adj": dge_2_pval_adj,
							"DGE_logfc": dge_2_logfc,
							"RF_A": rf_obj.accuracy,
							"RF_AUC": rf_obj.auc_score,
							"RF_BP": rf_obj.rf_params,
						}

						images[str(primary)+'-'+str(secondary)] = rf_obj.clustering_details

					# only run if there are more than 2 categories
					if len(self.categories) > 2:

						#get the unique predicted_tertiary_costa for the secondary
						tertiaries = self.df[(self.df[self.categories[0]] == primary) & (self.df[self.categories[1]] == secondary)][self.categories[2]].unique()

						#iterate through the tertiaries
						for tertiary in tertiaries:

							try:

								# change the focus of the rf_obj
								rf_obj.change_focus(primary=primary,secondary=secondary,tertiary=tertiary,serotype_of_interest=gene_of_interest,min_cells=self.minCells,min_serotype_cells=self.minSeroTypeCells)

								# dge data dictionary
								dge_3 = rf_obj.de_df.set_index('1_name')['0_name'].to_dict()
								dge_3_pval_adj = rf_obj.de_df[['1_pval_adj','0_pval_adj']].values.tolist()
								dge_3_logfc = rf_obj.de_df[['1_logfc','0_logfc']].values.tolist()

								# rf dictionary
								rf_3 = rf_obj.fi_df[:50].set_index('Feature')['Importance'].to_dict()

							except Exception as e:
								warnings.warn(f"Warning (4): {node} {primary}:{secondary}:{tertiary}:{gene_of_interest} {e}")
								dge_3 = {"Error": primary+":"+secondary+":"+tertiary}
								rf_3 = {"Error": primary+":"+secondary+":"+tertiary}
								dge_3_pval_adj = {"Error": primary+":"+secondary+":"+tertiary}
								dge_3_logfc = {"Error": primary+":"+secondary+":"+tertiary}
								rf_obj.clustering_details = ""


							# create the tertiary category key
							if node:

								data_structure[node][primary][secondary][tertiary] = {
									"RF": rf_3,
									"DGE": dge_3,
									"DGE_pval_adj": dge_3_pval_adj,
									"DGE_logfc": dge_3_logfc,
									"RF_A": rf_obj.accuracy,
									"RF_AUC": rf_obj.auc_score,
									"RF_BP": rf_obj.rf_params,
								}
								
								images[str(node)+'-'+str(primary)+'-'+str(secondary)+'-'+str(tertiary)] = rf_obj.clustering_details

							else:

								data_structure[primary][secondary][tertiary] = {
									"RF": rf_3,
									"DGE": dge_3,
									"DGE_pval_adj": dge_3_pval_adj,
									"DGE_logfc": dge_3_logfc,
									"RF_A": rf_obj.accuracy,
									"RF_AUC": rf_obj.auc_score,
									"RF_BP": rf_obj.rf_params,
								}

								images[str(primary)+'-'+str(secondary)+'-'+str(tertiary)] = rf_obj.clustering_details


		return data_structure



	def safe_gt(self, val, threshold):
		try:
			return float(val) > threshold
		except (ValueError, TypeError):
			return False

	def add_key_to_dict(self, d, key, value, default_value=1):
		if key not in d:
			d[key] = {k: default_value for k in value}

		return d

	def is_pos(self,v):
		try:
			return float(v) > 0
		except (ValueError, TypeError):
			return False


	def find_pepper_results(self):

		ignore_nodes = ['RF', 'DGE', 'DGE_pval_adj', 'DGE_logfc', 'RF_A', 'RF_AUC', 'RF_BP']

		pepper_results = {}

		#iterate all nodes of the tree
		for node in self.data_structure:

			#print(f"===============\n{node} \n===============")

			#COLUMNS
			#are there any other nodes at this level with the same name?
			duplicates_cols1 = [n for n in self.data_structure if n == node]
			if len(duplicates_cols1) > 1:
				#print(f"\tDuplicates: {duplicates_cols1}")
				
				#add to the pepper_results dict
				self.add_key_to_dict(pepper_results, node, duplicates_cols1)


			#ROWS
			#are there any nodes at this level with the same name minus the last character?
			duplicates_rows = [n for n in self.data_structure if n[:-1] == node[:-1]]
			if len(duplicates_rows) > 1:

				#check the node.RF of the duplicates to see if any of the keys are the same and values  > 0 
				#duplicates_rf = [list(filter(lambda x: self.data_structure[n]["RF"][x] > 0, self.data_structure[n]["RF"].keys())) for n in duplicates_rows]		
				
				duplicates_rf = [
					[k for k, v in self.data_structure[n]["RF"].items() if self.is_pos(v)]
					for n in duplicates_rows
				]
				duplicates_rf = set.intersection(*map(set, duplicates_rf))
				
				#check the node.DGE of the duplicate, but just look for the keys. They don't need to be > 0
				duplicates_dge = [list(self.data_structure[n]["DGE"].keys()) for n in duplicates_rows]
				duplicates_dge = set.intersection(*map(set, duplicates_dge))
				
				#which appear in both RF and DGE
				duplicates_both = [x for x in duplicates_rf if x in duplicates_dge]

				#print(f"\tRows: {duplicates_rows}: peppers1: {duplicates_both}")
				
				#add to the pepper_results dict
				self.add_key_to_dict(pepper_results, node, duplicates_both)
			

			#show all the keys this node that are in self.data_structure.keys()
			keys = [k for k in self.data_structure.keys() if k in self.data_structure[node].keys()]

			#iterate all of the node keys
			if len(duplicates_rows) > 0:
				for k in keys:
					list_buffer_rf = []
					list_buffer_dge = []
					for row in duplicates_rows:
						
						list_buffer_rf.append([
							x for x in self.data_structure[row][k]['RF']
							if self.safe_gt(self.data_structure[row][k]['RF'][x], 0)
						])

						#same for DGE, but doesn't need to be there
						list_buffer_dge.append([
							x for x in self.data_structure[row][k]['DGE']
						])
						
					
					#intersection of all values in each list
					list_buffer_rf = list(set.intersection(*map(set, list_buffer_rf)))
					list_buffer_dge = list(set.intersection(*map(set, list_buffer_dge)))
					
					#which appear in both RF and DGE
					list_buffer_both = [x for x in list_buffer_rf if x in list_buffer_dge]

					#print(f"\n\t{k}:{duplicates_rows}: peppers2: {list_buffer_both}")
					
					#add to the pepper_results dict
					self.add_key_to_dict(pepper_results, node+"-"+k, list_buffer_both)
					
					duplicates_cols2 = [n for n in self.data_structure[node].keys() if n[:-1] == k[:-1]]
					for col in duplicates_cols2:
						for key2 in self.data_structure[node][col].keys():
							if key2 not in ignore_nodes:
								list_buffer_rf_2 = []
								list_buffer_dge_2 = []

								for col2 in duplicates_cols2:
									try:
										list_buffer_rf_2.append([
											x for x in self.data_structure[node][col2][key2]['RF']
											if self.safe_gt(self.data_structure[node][col2][key2]['RF'][x], 0)
										])

										list_buffer_dge_2.append([
											x for x in self.data_structure[node][col2][key2]['DGE']
										])
									except KeyError:
										continue

								#intersection of all values in each list
								list_buffer_rf_2 = list(set.intersection(*map(set, list_buffer_rf_2)))
								list_buffer_dge_2 = list(set.intersection(*map(set, list_buffer_dge_2)))

								#which appear in both RF and DGE
								list_buffer_both_2 = [x for x in list_buffer_rf_2 if x in list_buffer_dge_2]
								
								#print(f"\t{node}:{col}:{key2}: {duplicates_cols2}: peppers3: {list_buffer_both_2}")
								
								#add to the pepper_results dict
								self.add_key_to_dict(pepper_results, node+"-"+col+"-"+key2, list_buffer_both_2)

		return pepper_results

	def get_unique_peppers(self, pepper_results):
		pepper_list = []
		for key in pepper_results.keys():
			for pepper in pepper_results[key]:
				if pepper not in pepper_list:
					pepper_list.append(pepper)
		return pepper_list

	def get_pepper_gene_details(self, pepper_list):
		pepper_genes = []
		mg = get_client('gene')
		for pepper in pepper_list:
			gene = mg.query(pepper, size=5)
			if gene["hits"]:
				gene_id = gene["hits"][0]["_id"]
				gene_details = str(mg.getgene(gene_id, fields=['name', 'symbol', 'summary','go'])).lower()
				if any(term in gene_details for term in ['plasmamembrane', 'plasma membrane', 'transmembrane', 'adhesion', 'receptor', 'extracellular']):
					pepper_genes.append(pepper)
					
		return pepper_genes

	def increase_pepper_heat(self, pepper_results, peppers):
		for key, value in pepper_results.items():
			for k, v in value.items():
				if k in peppers:
					pepper_results[key][k] = pepper_results[key][k] + 1

		return pepper_results

