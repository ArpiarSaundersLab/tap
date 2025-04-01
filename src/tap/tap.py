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
import importlib.resources as pkg_resources
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
	
		# self.loadData()
		# self.runTimeEstimate()

		# if self.runCellTypist == True:
		# 	self.generateCellTypistPredictions()

		# self.generateHeatmapMetadata()

		# if self.excludeMarkers == True:
		# 	self.determineMarkers()

		# if self.mapOnly == True:
		# 	self.generateCSV()
		# else:
		# 	self.generateCSV()
		# 	self.generateJsonFile()

		# if minify == True:	
		# 	self.minifyHtml()


	def copyTapTemplate(self):		
		with pkg_resources.path('tap.templates', 'tap.html') as template_path:
			shutil.copy(template_path, self.outputPath + self.outputName)
