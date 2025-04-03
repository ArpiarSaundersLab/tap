import scanpy as sc
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import tables
import copy
import warnings
import gc
import json
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
warnings.filterwarnings("ignore")
sns.set(style="whitegrid")
sc.settings.verbosity = 0

class FeatureSelection:

	def __init__(self, molotAAV_object, serotype_of_interest, primary=None, exclude_category_names=[], 
				secondary=None, tertiary=None, min_cells=100, min_serotype_cells=100, 
				permute_feature_importance=False, permutation_repeats=2, threads=8,
				cluster_method='threshold', cluster_threshold=1, showDE=True, showRF=True, serotype_list=[],
				category_observations=[], showPlots=True, exclude_genes=[], deMethod="wilcoxon",
				rfHyperParameterTune=False,rfHyperParameterIterations=50, rfType="classifier", outputPath=None, balance=None,
				excludeMarkers=False,markers=[],removeOutliers=True,totalIterations=0, currentIteration=1,
				clusterThresholdGreaterThanOrEqual=None, clusterThresholdLessThanOrEqual=None):
		self.molotAAV_object = molotAAV_object	
		self.molotAAV_object.adata = molotAAV_object
		self.serotype_of_interest = serotype_of_interest.lower()
		self.primary = primary
		self.secondary = secondary
		self.tertiary = tertiary
		self.categories = [self.primary, self.secondary, self.tertiary]
		self.min_cells = min_cells
		self.min_serotype_cells = min_serotype_cells
		self.threads = threads
		self.cluster_method = cluster_method
		self.cluster_threshold = cluster_threshold
		self.clusterThresholdGreaterThanOrEqual = clusterThresholdGreaterThanOrEqual
		self.clusterThresholdLessThanOrEqual = clusterThresholdLessThanOrEqual
		self.showDE = showDE
		self.showRF = showRF
		self.deMethod = deMethod
		self.showPlots = showPlots
		self.outputPath = outputPath
		self.exclude_genes = exclude_genes
		self.category_observations = category_observations
		self.exclude_category_names = exclude_category_names
		self.rfHyperParameterTune = rfHyperParameterTune
		self.rfHyperParameterIterations = rfHyperParameterIterations
		self.rfType = rfType
		self.balance = balance
		self.accuracy = "0%"
		self.auc_score = 0
		self.rf_params= ""
		self.excludeMarkers = excludeMarkers
		self.markers = markers
		self.removeOutliers = removeOutliers
		self.totalIterations = totalIterations
		self.currentIteration = currentIteration
		self.clustering_details = {}
		self.fi_df = pd.DataFrame(columns=['Feature', 'Importance'])
		self.fip_df = pd.DataFrame(columns=['Feature', 'Importance'])
		self.de_df = pd.DataFrame()
		self.molotAAV_object.serotype_list = serotype_list
		self.molotAAV_object_copy = copy.deepcopy(self.molotAAV_object)
		self.convert_to_upper()
		if self.removeOutliers == True:
			self.preprocess()
		self.add_log10()
		self.filter_by_categories()
		self.print_run_details()

		if self.removeOutliers == True:
			self.remove_outliers()
		else:
			self.molotAAV_object_filtered = self.molotAAV_object_processed.adata
		
		if self.check_shape_count() == False:
			print(f"Too few cells in the cell type to proceed.")
			return
		
		if self.check_shape_serotype() == False:
			print(f"Not enough cells with the {self.serotype_of_interest} to proceed.")
			return

		if self.check_labels() == False:
			print(f"Not enough unique labels in the {self.serotype_of_interest} to proceed.")
			return

		if self.cluster_method == 'gaussian':
			self.run_gmm()
		else:
			self.cluster_by_threshold()

		if self.showRF == True:
			self.run_rf(permute_feature_importance=permute_feature_importance, permutation_repeats=permutation_repeats)
			self.print_feature_importances()

		self.differential_expression()
		self.plot_serotype_vs_umi_training()

	def change_focus(self, serotype_of_interest, primary=None, secondary=None, tertiary=None, min_cells=100, min_serotype_cells=100, permute_feature_importance=False, permutation_repeats=2):
		"""
		Change the focus of the analysis to a different serotype, primary, secondary, or tertiary.
		"""
		self.clear_memory()
		self.serotype_of_interest = serotype_of_interest.lower()
		self.primary = primary
		self.secondary = secondary
		self.tertiary = tertiary
		self.categories = [self.primary, self.secondary, self.tertiary]
		self.min_cells = min_cells
		self.min_serotype_cells = min_serotype_cells
		self.accuracy = "0%"
		self.auc_score = 0
		self.add_log10()
		self.filter_by_categories()
		self.currentIteration += 1
		self.print_run_details()
		
		
		if self.removeOutliers == True:
			self.remove_outliers()
		else:
			self.molotAAV_object_filtered = self.molotAAV_object_processed.adata

		if self.check_shape_count() == False:
			print(f"Too few cells in the cell type to proceed.")
			return
		
		if self.check_shape_serotype() == False:
			print(f"Not enough cells with the {self.serotype_of_interest} to proceed.")
			return

		if self.check_labels() == False:
			print(f"Not enough unique labels in the {self.serotype_of_interest} to proceed.")
			return

		if self.cluster_method == 'gaussian':
			self.run_gmm()
		else:
			self.cluster_by_threshold()
		
		if self.showRF == True:
			self.run_rf(permute_feature_importance=permute_feature_importance, permutation_repeats=permutation_repeats)
			self.print_feature_importances()

		self.differential_expression()
		self.plot_serotype_vs_umi_training()
		
	def print_run_details(self):
		print("("+str(self.currentIteration)+"/"+str(self.totalIterations)+") "+str(self.serotype_of_interest) + ":"+str(self.primary) + ":"+str(self.secondary) + ":"+str(self.tertiary))

	def convert_to_upper(self):	
		self.molotAAV_object_copy.adata.var_names = [x.upper() for x in self.molotAAV_object_copy.adata.var_names]
		self.molotAAV_object_copy.adata.var_names_make_unique()

	def preprocess(self):
		if self.exclude_category_names != []:
			self.molotAAV_object_copy.adata = self.molotAAV_object_copy.adata[~self.molotAAV_object_copy.adata.obs[self.primary].isin(self.exclude_category_names), :]
		self.molotAAV_object_copy.adata = self.molotAAV_object_copy.adata[self.molotAAV_object_copy.adata.obs["total_counts"] < self.molotAAV_object_copy.adata.obs["total_counts"].quantile(0.98), :]
		self.molotAAV_object_copy.adata = self.molotAAV_object_copy.adata[self.molotAAV_object_copy.adata.obs["infection_count"] < self.molotAAV_object_copy.adata.obs["infection_count"].quantile(0.98), :]

	def add_log10(self):
		self.molotAAV_object_copy.adata.obs["infection_count_log10"] = np.log10(self.molotAAV_object_copy.adata.obs["infection_count"])
		self.molotAAV_object_copy.adata.obs["infection_count_log10"][np.isinf(self.molotAAV_object_copy.adata.obs["infection_count_log10"])] = 0
		self.molotAAV_object_copy.adata.obs["total_counts_log10"] =  np.log10(self.molotAAV_object_copy.adata.obs["total_counts"])
		self.molotAAV_object_copy.adata.obs["total_counts_log10"][np.isinf(self.molotAAV_object_copy.adata.obs["total_counts_log10"])] = 0

		for serotype in self.molotAAV_object_copy.serotype_list:
			serotype = serotype.upper()
			# does the serotype exist in the adata object?
			if serotype in self.molotAAV_object_copy.adata.var_names:
				self.molotAAV_object_copy.adata[self.molotAAV_object_copy.adata[:, [serotype]].X == 0][:,serotype].X = 0.1
				self.molotAAV_object_copy.adata.obs[str(serotype).lower() + "_log10"] = np.log10(self.molotAAV_object_copy.adata[:, [serotype]].to_df())
				self.molotAAV_object_copy.adata.obs[str(serotype).lower()] =self.molotAAV_object_copy.adata[:, [serotype]].to_df()
 
	def filter_by_categories(self):
		self.molotAAV_object_processed = copy.deepcopy(self.molotAAV_object_copy)
		inc = 0
		for category in self.category_observations:
			if self.categories[inc] != None:
				self.molotAAV_object_processed.adata = self.molotAAV_object_processed.adata[self.molotAAV_object_processed.adata.obs[category]==self.categories[inc], :]
			inc += 1
	
	def check_labels(self):

		if self.molotAAV_object_processed.adata.obs[self.serotype_of_interest].nunique() < 2:
			print(f"There must be at least 2 unique labels in the {self.serotype_of_interest} column.")
			return False
		else:
			return True

	def check_shape_count(self):
		# if the shape of the adata object is less than min_cells, then return false and stop the analysis.
		if self.molotAAV_object_processed.adata.shape[0] < self.min_cells:
			print("Only "+str(self.molotAAV_object_processed.adata.shape[0])+" cells found with the cell type of interest, "+str(self.primary))
			return False
		else:
			return True
	
	def check_shape_serotype(self):	
		# does the serotype of interest have log10 in the name?
		if "_log10" in self.serotype_of_interest:
			serotype_count = self.molotAAV_object_processed.adata.obs[self.serotype_of_interest.lower()][self.molotAAV_object_processed.adata.obs[self.serotype_of_interest.lower()] > 0.1].shape[0]
		else:
			serotype_count = self.molotAAV_object_processed.adata.obs[self.serotype_of_interest.lower()][self.molotAAV_object_processed.adata.obs[self.serotype_of_interest.lower()] > 0].shape[0]
		# if there are less than  min_serotype_cells with the serotype of interest, then return false and stop the analysis.
		if serotype_count < self.min_serotype_cells:
			print("Only "+str(serotype_count)+" cells with the serotype of interest, "+str(self.serotype_of_interest)+" found.")
			return False
		else:
			return True
		
	def remove_outliers(self, lower_percentile=0, upper_percentile=99):
		lower_percentile = np.percentile(self.molotAAV_object_processed.adata.obs[self.serotype_of_interest], lower_percentile)
		upper_percentile = np.percentile(self.molotAAV_object_processed.adata.obs[self.serotype_of_interest], upper_percentile)
		if np.isnan(lower_percentile):
			lower_percentile = 0
		if np.isnan(upper_percentile):
			upper_percentile = 100
		self.molotAAV_object_filtered = self.molotAAV_object_processed.adata[(self.molotAAV_object_processed.adata.obs[self.serotype_of_interest] >= lower_percentile) & (self.molotAAV_object_processed.adata.obs[self.serotype_of_interest] < upper_percentile)]

	def generateFigureTitle(self):
		name =  str(self.serotype_of_interest)
		if self.primary != None:
			name += " | "+str(self.primary)
		if self.secondary != None:
			name += " | "+str(self.secondary)
		if self.tertiary != None:
			name += " | "+str(self.tertiary)

		if name == "infection_count" or name == "infection_count_log10":
			name = "All Infection Counts"

		return str(name)

	def generateFigureFilename(self):
		filename = ""
		if self.serotype_of_interest != None:
			if self.serotype_of_interest != "infection_count":
				filename += str(self.serotype_of_interest)
		if self.primary != None:
			if filename == "":
				filename += str(self.primary)
			else:
				filename += "-"+str(self.primary)
		if self.secondary != None and self.serotype_of_interest != None:
			filename += "-"+str(self.secondary)
		if self.tertiary != None and self.secondary != None and self.serotype_of_interest != None:
			filename += "-"+str(self.tertiary)

		if filename == "":
			filename = "All"

		return str(filename).lower()

	def plot_serotype_vs_umi(self):
		plt.figure(figsize=(4, 4))
		g=sns.scatterplot(x=self.molotAAV_object_processed.adata.obs["total_counts_log10"], 
						y=self.molotAAV_object_processed.adata.obs[self.serotype_of_interest])
		g.set(title="All Cells\n"+self.generateFigureTitle())
		plt.tight_layout()
		plt.savefig(self.outputPath+'/images/'+self.generateFigureFilename()+'-library.png', dpi=80, facecolor='w', edgecolor='w')
		#plt.show()
		plt.clf()

	def plot_serotype_vs_umi_filtered(self, plot=False):
		plt.figure(figsize=(4, 4))
		g=sns.scatterplot(x=self.molotAAV_object_filtered.obs["total_counts_log10"], 
						y=self.molotAAV_object_filtered.obs[self.serotype_of_interest])
		g.set(title="Filtered Cells\n"+self.generateFigureTitle())
		plt.tight_layout()
		plt.savefig(self.outputPath+'/images/'+self.generateFigureFilename()+'-filtered.png', dpi=80, facecolor='w', edgecolor='w')
		#plt.show()
		plt.clf()

	def plot_serotype_vs_umi_training(self):
		fig = plt.figure(figsize=(4, 4))
		g = sns.scatterplot(x=self.molotAAV_object_filtered_training.obs["total_counts_log10"],
							y=self.molotAAV_object_filtered_training.obs[self.serotype_of_interest],
							hue=self.molotAAV_object_filtered_training.obs['infection_status'])
		g.set(title="Clustered\n"+self.generateFigureTitle())
		g.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6})
		plt.tight_layout()
		#plt.savefig(self.outputPath+'/images/'+self.generateFigureFilename()+'-clusters.png', dpi=80, facecolor='w', edgecolor='w')
		plt.savefig(self.outputPath+'clusters.png', dpi=80, facecolor='w', edgecolor='w')
		plt.clf()
		plt.close()

		import base64
		import os

		#save as base64 to use in the web application
		with open(self.outputPath+'/clusters.png', "rb") as image_file:
			encoded_string = base64.b64encode(image_file.read())
			encoded_string = encoded_string.decode("utf-8")
			self.clustering_details = encoded_string

		#delete the file
		os.remove(self.outputPath+'/clusters.png')


	def generate_cluster_details(self):
		self.clustering_details = [self.molotAAV_object_filtered_training.obs["total_counts_log10"].to_list(),
							self.molotAAV_object_filtered_training.obs[self.serotype_of_interest].to_list(),
							self.molotAAV_object_filtered_training.obs['infection_status'].to_list()]

		

	def run_gmm(self,n_components=2,max_iter=100,covariance_type="spherical"):
		#extract the features to cluster
		X = self.molotAAV_object_filtered.obs[["total_counts_log10", self.serotype_of_interest]]
		
		if X.empty:
			print(self.molotAAV_object_filtered.obs)
			raise ValueError("The feature set X is empty. Please check your data loading and preprocessing steps.")

		gmm = GaussianMixture(n_components=n_components,max_iter=max_iter,covariance_type=covariance_type).fit(X)
		
		#make predictions of each cluster / infection status
		self.molotAAV_object_filtered.obs['infection_status'] = gmm.predict(X)

		#get the probabilities of each instance in X
		probs = gmm.predict_proba(X)

		#store the probabilities in the DataFrame
		self.molotAAV_object_filtered.obs['prob_cluster_0'] = probs[:, 0]
		self.molotAAV_object_filtered.obs['prob_cluster_1'] = probs[:, 1]

		#select only cells with cluster_1 or cluster_0 greater than 95% and add them to a new object
		self.molotAAV_object_filtered_training = self.molotAAV_object_filtered[(self.molotAAV_object_filtered.obs['prob_cluster_0'] > 0.95) | (self.molotAAV_object_filtered.obs['prob_cluster_1'] > 0.95)]


	def cluster_by_threshold(self):
		self.molotAAV_object_filtered_training = self.molotAAV_object_filtered.copy()
		obs = self.molotAAV_object_filtered_training.obs

		if self.clusterThresholdGreaterThanOrEqual is not None and self.clusterThresholdLessThanOrEqual is not None:
			#create a new column for infection_status
			infection_status = pd.Series(index=obs.index, dtype='float')

			#assign values based on the conditions
			infection_status[obs[self.serotype_of_interest] >= self.clusterThresholdGreaterThanOrEqual] = 1
			infection_status[obs[self.serotype_of_interest] <= self.clusterThresholdLessThanOrEqual] = 0
			infection_status[
				(obs[self.serotype_of_interest] > self.clusterThresholdLessThanOrEqual) & 
				(obs[self.serotype_of_interest] < self.clusterThresholdGreaterThanOrEqual)
			] = None  #assign None for the ambiguous range

			#assign back to obs
			self.molotAAV_object_filtered_training.obs['infection_status'] = infection_status

			#remove rows with None if required
			self.molotAAV_object_filtered_training = self.molotAAV_object_filtered_training[
				~self.molotAAV_object_filtered_training.obs['infection_status'].isna()
			]
		else:
			#single threshold for clustering
			self.molotAAV_object_filtered_training.obs['infection_status'] = np.where(
				self.molotAAV_object_filtered_training.obs[self.serotype_of_interest] >= self.cluster_threshold,
				1,0
			)

	def print_infectivity_counts(self):
		print(self.molotAAV_object_filtered_training.obs['infection_status'].value_counts())

	def run_rf(self, n_estimators=20, max_depth=2, n_jobs=8, max_features=None, min_samples_split=2, permute_feature_importance=False, permutation_repeats=2):
		#filter out the serotype of interest from the training data
		self.molotAAV_object_filtered_training = self.molotAAV_object_filtered_training[:, ~self.molotAAV_object_filtered_training.var_names.isin([x.upper() for x in self.molotAAV_object_copy.serotype_list])]
		self.molotAAV_object_filtered_training = self.molotAAV_object_filtered_training[:, ~self.molotAAV_object_filtered_training.var_names.isin([x.upper() for x in self.exclude_genes])]

		#we need to filter out any marker genes from the list
		if self.excludeMarkers == True:
			if self.primary != None and self.secondary != None and self.tertiary != None:
				exclude_markers = self.markers[self.tertiary]
			elif self.primary != None and self.secondary != None:
				exclude_markers = self.markers[self.secondary]
			elif self.primary != None:
				exclude_markers = self.markers[self.primary]
			else:
				exclude_markers = {}
				for key in self.markers:				
					if "names" not in exclude_markers:
						exclude_markers["names"] = []
					exclude_markers["names"].append(self.markers[key]["names"])
				exclude_markers["names"] = [item for sublist in exclude_markers["names"] for item in sublist]

			self.molotAAV_object_filtered_training = self.molotAAV_object_filtered_training[:, ~self.molotAAV_object_filtered_training.var_names.isin([x.upper() for x in exclude_markers["names"]])]

		#add obs that's a category of infected or not infected
		self.molotAAV_object_filtered_training.obs["infection_status_cat"] = self.molotAAV_object_filtered_training.obs["infection_status"].astype('category')		

		#set the x training and y labels
		X_train = self.molotAAV_object_filtered_training.X
		y_train = self.molotAAV_object_filtered_training.obs["infection_status"].values		

		#both classes must be present in the training data
		if len(np.unique(y_train)) < 2:
			return

		# balance the data according to the balance parameter
		if self.balance is not None:
			if self.balance.lower() == "smote":
				from imblearn.over_sampling import SMOTE
				smote = SMOTE(random_state=22)
				X_resampled, y_resampled = smote.fit_resample(X_train, y_train)
				X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=22)
			elif self.balance.lower() == "randomoversampler":
				from imblearn.over_sampling import RandomOverSampler
				ros = RandomOverSampler(random_state=22)
				X_resampled, y_resampled = ros.fit_resample(X_train, y_train)
				X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=22)
			elif self.balance.lower() == "adasyn":
				from imblearn.over_sampling import ADASYN
				ada = ADASYN(random_state=22)
				X_resampled, y_resampled = ada.fit_resample(X_train, y_train)
				X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=22)
			elif self.balance.lower() == "borderlinesmote":
				from imblearn.over_sampling import BorderlineSMOTE
				bsmote = BorderlineSMOTE(random_state=22)
				X_resampled, y_resampled = bsmote.fit_resample(X_train, y_train)
				X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=22)
			elif self.balance.lower() == "kmeanssmote":
				from imblearn.over_sampling import KMeansSMOTE
				ksmote = KMeansSMOTE(random_state=22)
				X_resampled, y_resampled = ksmote.fit_resample(X_train, y_train)
				X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=22)
			elif self.balance.lower() == "svmsmote":
				from imblearn.over_sampling import SVMSMOTE
				ssmote = SVMSMOTE(random_state=22)
				X_resampled, y_resampled = ssmote.fit_resample(X_train, y_train)
				X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=22)
			else:
				print("Invalid balance parameter. Please use one of the following: smote, randomoversampler, adasyn, borderlinesmote, kmeanssmote, svmsmote.")
				return

		else:
			#split the training data into training (x%), testing (x%) sets
			X_train, X_test, y_train, y_test = train_test_split(X_train, y_train, test_size=0.2, random_state=22)

		#use a random forest classifer or regressor
		if self.rfType == "classifier":

			if self.rfHyperParameterTune == False:
				#initialize a random forest classifier
				self.rfc = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, 
											n_jobs=self.threads, max_features=max_features, min_samples_split=min_samples_split, random_state=22)

				#fit the classifier to the training data
				self.rfc.fit(X_train, y_train)

				#make predictions on the testing data
				y_pred = self.rfc.predict(X_test)
				self.rf_params = {"n_estimators":n_estimators, "max_depth":max_depth, "n_jobs":self.threads, "max_features":max_features, "min_samples_split":min_samples_split}
				self.rf_params = json.dumps(self.rf_params)

				#evaluate the accuracy of the model on the testing data
				accuracy = accuracy_score(y_test, y_pred)
				self.accuracy = str(round(accuracy*100,2))+"%"
				self.auc_score = roc_auc_score(y_test, y_pred)
				
			else:
				n_estimators = [10,20,40,60,80,100]
				max_features = ['auto', 'sqrt', None]
				max_depth = [2,3,4,5,6,7,8,9,10,11,12,None]
				min_samples_split = [2,3,4, 5,7,8,9,10]
				min_samples_leaf = [1, 2, 3, 4, 5]

				#create the grid of parameters
				random_grid = {'n_estimators': n_estimators,
					'max_features': max_features,
					'max_depth': max_depth,
					'min_samples_split': min_samples_split,
					'min_samples_leaf': min_samples_leaf
					}

				#use the random grid to search for best hyperparameters
				#rf = RandomForestRegressor()
				rf = RandomForestClassifier(random_state=22)
				rf_random = RandomizedSearchCV(estimator=rf, param_distributions=random_grid, 
											n_iter=self.rfHyperParameterIterations, 
											cv=3, 
											verbose=0, 
											random_state=22, 
											n_jobs=self.threads)
				rf_random.fit(X_train, y_train)
				
				self.rf_params = rf_random.best_params_
				self.rf_params = json.dumps(self.rf_params)
				self.rfc = rf_random.best_estimator_
				y_pred = self.rfc.predict(X_test)
				y_pred_binary = (y_pred > 0.5).astype(int)
				accuracy = accuracy_score(y_test, y_pred_binary)
				self.accuracy = str(round(accuracy*100,2))+"%"
				self.auc_score = roc_auc_score(y_test, y_pred)
				
		elif self.rfType == "regressor":

			if self.rfHyperParameterTune == False:
				#initialize a random forest regressor
				self.rfc = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth, 
											n_jobs=self.threads, max_features=max_features, min_samples_split=min_samples_split, random_state=22)

				#fit the regressor to the training data
				self.rfc.fit(X_train, y_train)

				#make predictions on the testing data
				y_pred = self.rfc.predict(X_test)
				self.rf_params = {"n_estimators":n_estimators, "max_depth":max_depth, "n_jobs":self.threads, "max_features":max_features, "min_samples_split":min_samples_split}
				self.rf_params = json.dumps(self.rf_params)

				#evaluate the accuracy of the model on the testing data
				# accuracy = accuracy_score(y_test, y_pred)
				# self.accuracy = f"{round(accuracy*100,2)}%"
				# self.auc_score = roc_auc_score(y_test, y_pred)

				y_pred_binary = (y_pred > 0.5).astype(int)
				accuracy = accuracy_score(y_test, y_pred_binary)
				self.accuracy = str(round(accuracy*100,2))+"%"
				self.auc_score = roc_auc_score(y_test, y_pred)


			else:
				n_estimators = [10,20,40,60,80,100]
				max_features = ['auto', 'sqrt', None]
				max_depth = [2,3,4,5,6,7,8,9,10,11,12,None]
				min_samples_split = [2,3,4, 5,7,8,9,10]
				min_samples_leaf = [1, 2, 3, 4, 5]

				#create the grid of parameters
				random_grid = {'n_estimators': n_estimators,
					'max_features': max_features,
					'max_depth': max_depth,
					'min_samples_split': min_samples_split,
					'min_samples_leaf': min_samples_leaf
					}

				#use the random grid to search for best hyperparameters
				rf = RandomForestRegressor(random_state=22)
				rf_random = RandomizedSearchCV(estimator=rf, param_distributions=random_grid,
											n_iter=self.rfHyperParameterIterations,
											cv=3,
											verbose=0,
											random_state=22,
											n_jobs=self.threads)
				rf_random.fit(X_train, y_train)

				self.rf_params = rf_random.best_params_
				self.rf_params = json.dumps(self.rf_params)
				self.rfc = rf_random.best_estimator_
				y_pred = self.rfc.predict(X_test)
				y_pred_binary = (y_pred > 0.5).astype(int)
				accuracy = accuracy_score(y_test, y_pred_binary)
				self.accuracy = str(round(accuracy*100,2))+"%"
				self.auc_score = roc_auc_score(y_test, y_pred)

		#permute feature importance
		if permute_feature_importance == True:
			#convert X_test and y_test to a numpy array if necessary
			if isinstance(X_test, sp.spmatrix):
				print("The matrix is sparse. Converting...")
				X_test = X_test.toarray()

			result = permutation_importance(
				self.rfc, X_test, y_test, n_repeats=permutation_repeats, random_state=22, n_jobs=self.threads
			)
			feature_names = self.molotAAV_object_filtered_training.var_names[self.molotAAV_object_filtered_training.var_names != "infection_status"]
			self.fip_df = pd.DataFrame({"Feature": feature_names, "Importance": result.importances_mean})
			self.fip_df = self.fip_df.sort_values(by="Importance", ascending=False)


	def print_feature_importances(self):

		#get the feature importances from the random forest model
		importances = self.rfc.feature_importances_

		#limit scores to 4 decimal places
		importances = [round(x,4) for x in importances]

		#get the names of the features
		feature_names = self.molotAAV_object_filtered_training.var_names[self.molotAAV_object_filtered_training.var_names != "infection_status"]

		#create a pandas dataframe to store the feature importances
		df = pd.DataFrame({"Feature": feature_names, "Importance": importances})

		#sort the dataframe by importance in descending order
		df = df.sort_values(by="Importance", ascending=False)
		
		#assign the dataframe to the class property
		self.fi_df = df


	def plot_serotype_vs_umi_predictions(self, plot=False):
		
		if plot == True:
			#subset the validation data from the adata object
			infection_validation = copy.deepcopy(self.molotAAV_object_copy.adata)
			infection_validation = infection_validation[:, ~infection_validation.var_names.isin([x.upper() for x in self.molotAAV_object_copy.serotype_list])]

			#extract the features from the validation data
			X_val = infection_validation[:, infection_validation.var_names != "infection_status"].X

			#make predictions on the validation data
			y_val_pred = self.rfc.predict(X_val)

			#add the predictions to the adata object
			infection_validation.obs["infection_status_pred"] = y_val_pred

			#view the plot of predictions
			g=sns.scatterplot(x=infection_validation.obs["total_counts_log10"], 
							y=infection_validation.obs[self.serotype_of_interest], 
							hue=infection_validation.obs['infection_status_pred'],
							s=30)
			g.set(xlabel='UMI Counts log(10)', ylabel=self.serotype_of_interest, title='Predictions')
			plt.show()
			plt.clf()

	def differential_expression(self, method='wilcoxon', category='infection_status', n_genes=1000, use_raw=False, plot=False):

		#filter out the serotype of interest from the training data
		self.molotAAV_object_filtered_training = self.molotAAV_object_filtered_training[:, ~self.molotAAV_object_filtered_training.var_names.isin([x.upper() for x in self.molotAAV_object_copy.serotype_list])]
		self.molotAAV_object_filtered_training = self.molotAAV_object_filtered_training[:, ~self.molotAAV_object_filtered_training.var_names.isin([x.upper() for x in self.exclude_genes])]

		#convert the infection status to a category
		self.molotAAV_object_filtered_training.obs[category] = self.molotAAV_object_filtered_training.obs[category].astype('category');
		
		#perform the differential expression test that was specified in deMethod
		sc.tl.rank_genes_groups(self.molotAAV_object_filtered_training, category, method=self.deMethod, use_raw=use_raw)

		names = pd.DataFrame(self.molotAAV_object_filtered_training.uns['rank_genes_groups']['names'])
		scores = pd.DataFrame(self.molotAAV_object_filtered_training.uns['rank_genes_groups']['scores'])
		logfoldchanges = pd.DataFrame(self.molotAAV_object_filtered_training.uns['rank_genes_groups']['logfoldchanges'])
		pvals = pd.DataFrame(self.molotAAV_object_filtered_training.uns['rank_genes_groups']['pvals'])
		pvals_adj = pd.DataFrame(self.molotAAV_object_filtered_training.uns['rank_genes_groups']['pvals_adj'])

		df_de = pd.concat([names.iloc[:n_genes], scores.iloc[:n_genes], logfoldchanges.iloc[:n_genes], pvals.iloc[:n_genes], pvals_adj.iloc[:n_genes]], axis=1)
		df_de.columns = ["0_name","1_name","0_score","1_score","0_logfc","1_logfc","0_pval","1_pval","0_pval_adj","1_pval_adj"]
		self.de_df = df_de

		#filter  de_df to pvals < 0.05
		self.de_df = self.de_df[(self.de_df["0_pval_adj"] <= 0.05) | (self.de_df["1_pval_adj"] <= 0.05)]

		#convert 0_pval_adj values in df to only 2 decimal places, but keep scientific notation
		self.de_df["0_pval_adj"] = self.de_df["0_pval_adj"].apply(lambda x: '%.2E' % x)
		self.de_df["1_pval_adj"] = self.de_df["1_pval_adj"].apply(lambda x: '%.2E' % x)
		
		#0_logfc to 2 decimal places
		self.de_df["0_logfc"] = self.de_df["0_logfc"].apply(lambda x: round(x,2))
		self.de_df["1_logfc"] = self.de_df["1_logfc"].apply(lambda x: round(x,2))


	def clear_memory(self):
		self.molotAAV_object_processed = None
		self.molotAAV_object_filtered = None
		self.molotAAV_object_filtered_training = None
		self.rfc = None
		gc.collect()
		self.fi_df = pd.DataFrame(columns=['Feature', 'Importance'])
		self.fip_df = pd.DataFrame(columns=['Feature', 'Importance'])
		self.de_df = pd.DataFrame(columns=["0_name","1_name","0_score","1_score","0_logfc","1_logfc","0_pval","1_pval","0_pval_adj","1_pval_adj"])
		#print("\n========================\nGarbage Cleanup\n========================\n")