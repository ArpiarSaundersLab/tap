import unittest
import tap as t
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

class TestBuildTAP(unittest.TestCase):
	def setUp(self):
		self.adata = sc.datasets.pbmc68k_reduced()

	def build_tap(self):
		parameters = {
			"filename" : "../assets/data/PCCM.h5ad",
			"name": "PCCM",
			"genes": ["AAV2","AAV9","Anc80"],
			"categories": ["cell_type"],
			"outputPath": ".", #current directory
			"outputName": "test.html",
			"runCellTypist" : False,
			"cellTypistModel": "BrainCellData_268_CellTypes.pkl",
			"excludeMarkers" : True,
			"mapOnly" : False,
			"showRF": True,
			"clusterMethod": "threshold",
			"clusterThreshold": 1,
		}

		results = t.TAP(**parameters)

		self.assertEqual(len(result), self.adata.shape[0])

if __name__ == '__main__':
	unittest.main()