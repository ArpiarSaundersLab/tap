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
import anndata
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
warnings.filterwarnings("ignore")
sns.set(style="whitegrid")
sc.settings.verbosity = 0

class TAP:

	colors = [(0.8, 0.8, 0.8)] + [(plt.cm.Reds(i / 255)) for i in range(80, 256)]
	cmap_custom = ListedColormap(colors)

	def __init__(self, name="Heatmap", filename=None, adataObject=None, categories=[], categoryNames=[], 
	genes=[], exclude=[], useRaw=False, outputPath=".", outputName="tap.html", cpus = 2, 
	mapOnly=False, showDE=True, showRF=False, showUMAP=True, showDetails=True, 
	clusterMethod="gaussian", clusterThreshold=0, useLog10=False, showPlots=False, 
	excludeGenes=[],deMethod="t-test", rfHyperParameterTune=False, rfHyperParameterIterations=50, 
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
		self.runTimeEstimate()
		if self.runCellTypist == True:
			self.generateCellTypistPredictions()
		self.generateHeatmapMetadata()
		if self.excludeMarkers == True:
			self.determineMarkers()
		if self.mapOnly == True:
			self.generateCSV()
		else:
			self.generateCSV()
			self.generateJsonFile()
		if minify == True:	
			self.minifyHtml()

	def copyTapTemplate(self):		
		shutil.copy(pkg_resources.resource_filename("tap", "templates")+"/tap.html", self.outputPath + self.outputName)

	def loadData(self):

		self.adata = sc.read(self.filename)

		#if user passes in a specific object name, use that object
		if self.adataObject:
			self.adata = self.adataObject

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

		#make var names unique
		self.adata.var_names_make_unique()

		#calculate the percent of counts from mitochondrial genes
		self.adata.var["mt"] = self.adata.var_names.str.startswith("mt-")
		sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
		self.adata = self.adata[self.adata.obs.pct_counts_mt < 10, :]

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
					self.totalIterations += len(self.adata[self.adata.obs[self.categories[0]] == category, :].obs[self.categories[1]].unique())
				if category_length > 2:
					for subcategory in self.adata[self.adata.obs[self.categories[0]] == category, :].obs[self.categories[1]].unique():
						self.totalIterations += len(self.adata[(self.adata.obs[self.categories[0]] == category) & (self.adata.obs[self.categories[1]] == subcategory), :].obs[self.categories[2]].unique())
		
		#calculate the total iterations based on the parameters
		self.totalIterations  = (self.totalIterations*len(self.genes))+len(self.genes)+1


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
		with open(self.outputPath+'metadata.txt', 'w') as file:
			file.write(json.dumps(metadata, indent=4))


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
					
					#keep only the markers with a score of 2.2 std deviations above the median.
					self.markers[catName] = self.markers[catName][self.markers[catName]["scores"] > self.markers[catName]["scores"].median() + 2.2*self.markers[catName]["scores"].std()]
					self.markers[catName] = self.markers[catName][self.markers[catName]["pvals_adj"] <= 0.05]
					self.markers[catName] = self.markers[catName].sort_values(by="scores", ascending=False)
					self.markers[catName] = self.markers[catName].reset_index(drop=True)
			except Exception as e:
				print(f"Error in determineMarkers2: {e}")
		
		del adata_copy
		gc.collect()


	def generateCSV(self, plot=False):

		#normalize by library size using Scanpy's normalize_total function
		sc.pp.normalize_total(self.adata, target_sum=1e4)
		
		#calculate the infection_count again after normalization
		self.adata.obs['infection_count'] = np.sum(self.adata[:, self.genes].to_df(), axis=1).astype(int)
		is_infected = self.adata.obs['infection_count'] > 0
		self.adata.obs['infected'] = is_infected.astype(str)

		#remove anything in the exclude list
		self.adata = self.adata[~self.adata.obs[self.categories[0]].isin(self.exclude), :]

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

		#normalize each of the serotypes the same way
		for serotype in self.genes:

			#extract the column for the serotype
			serotype_data = self.adata[:, serotype].X.toarray().flatten()

			#find the maximum value in the serotype data
			max_value = serotype_data.max().astype(float)
			
			#NEW The normalization value should be the max library size.
			normalization_value = 1e6
			
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


	