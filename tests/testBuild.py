import unittest
import tap as t
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
import os

class TestBuild(unittest.TestCase):
	def setUp(self):
		warnings.filterwarnings('ignore')
		self.adata = sc.AnnData(X=np.random.randint(0,100, size=(1500,250)), 
							obs=pd.DataFrame(index=[f"cell_{i}" for i in range(1500)]), 
							var=pd.DataFrame(data=[f"gene_{i}" for i in range(250)],columns=["gene_name"]))
		self.adata.var.index = [f"gene_{i}" for i in range(self.adata.shape[1])]
		self.adata.obs["cell_type"] = np.random.choice(["A","B"], self.adata.shape[0])
		#adds de gene. gene_1 should always have gene_4 as a diff expressed gene.
		self.adata[self.adata[:,"gene_1"].X >=50, 5].X = 0

		#simple preprocessing
		self.adata.raw = self.adata.copy() #raw layer required
		sc.pp.normalize_total(self.adata, target_sum=1e4)
		sc.pp.log1p(self.adata)
		sc.pp.pca(self.adata, n_comps=50)
		sc.pp.neighbors(self.adata, n_neighbors=10)
		sc.tl.umap(self.adata)
		sc.tl.leiden(self.adata, resolution=0.5)

	def test_build_tap_gaussian(self):
		print("Testing TAP with gaussian clustering...")
		parameters = {
			"adataObject" :self.adata,
			"name": "Unittest",
			"genes": ["gene_1","gene_2"],
			"categories": ["cell_type"],
			"outputPath": ".", 
			"outputName": "unittest.html",
			"clusterMethod": "gaussian",
		}
		results = t.TAP(**parameters)
		success_checks = 0
		try:
			if results.data_structure['gene_1']['A']['RF']['GENE_5']:
				success_checks += 1
		except KeyError:
			print("GENE_5 not found in RF")
		try:
			if results.data_structure['gene_1']['A']['DGE']['GENE_5']:
				success_checks += 1
		except KeyError:
			print("GENE_5 not found in DGE")
		try:
			if len(results.data_structure['gene_2']['A']['DGE']) == 0:
				success_checks += 1
		except KeyError:
			print("gene_2 DGE >=0 is incorrect")

		self.assertEqual(3, success_checks)

	def test_build_tap_threshold(self):
		print("Testing TAP with threshold clustering...")
		parameters = {
			"adataObject" :self.adata,
			"name": "unittest",
			"genes": ["gene_1","gene_2"],
			"categories": ["cell_type"],
			"outputPath": ".", 
			"outputName": "unittest.html",
			"removeOutliers": True,
			"clusterMethod": "threshold",
			"clusterThreshold": 50,
			"minify": False,
		}
		results = t.TAP(**parameters)
		success_checks = 0
		try:
			if results.data_structure['gene_1']['A']['RF']['GENE_5']:
				success_checks += 1
		except KeyError:
			print("GENE_5 not found in RF")
		try:
			if results.data_structure['gene_1']['A']['DGE']['GENE_5']:
				success_checks += 1
		except KeyError:
			print("GENE_5 not found in DGE")
		try:
			if len(results.data_structure['gene_2']['A']['DGE']) == 0:
				success_checks += 1
		except KeyError:
			print("gene_2 DGE >=0 is incorrect")

		self.assertEqual(3, success_checks)


	def test_build_tap_threshold_buffer(self):
		print("Testing TAP with threshold buffer clustering...")
		parameters = {
			"adataObject" :self.adata,
			"name": "unittest",
			"genes": ["gene_1","gene_2"],
			"categories": ["cell_type"],
			"outputPath": ".", 
			"outputName": "unittest.html",
			"removeOutliers": True,
			"clusterMethod": "threshold",
			#"clusterThreshold": 50,
			"clusterThresholdGreaterThanOrEqual": 50,
			"clusterThresholdLessThanOrEqual": 30,
			"minify": False,
		}
		results = t.TAP(**parameters)
		success_checks = 0
		try:
			if results.data_structure['gene_1']['A']['RF']['GENE_5']:
				success_checks += 1
		except KeyError:
			print("GENE_5 not found in RF")
		try:
			if results.data_structure['gene_1']['A']['DGE']['GENE_5']:
				success_checks += 1
		except KeyError:
			print("GENE_5 not found in DGE")
		try:
			if len(results.data_structure['gene_2']['A']['DGE']) == 0:
				success_checks += 1
		except KeyError:
			print("gene_2 DGE >=0 is incorrect")

		self.assertEqual(3, success_checks)

	def tearDown(self):
		if os.path.exists("unittest.html"):
			os.remove("unittest.html")
		if os.path.exists("metadata.txt"):
			os.remove("metadata.txt")

if __name__ == '__main__':
	unittest.main()