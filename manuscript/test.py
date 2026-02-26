import scanpy as sc
import celltypist
import pandas as pd
import numpy as np
import time

#opens the Cortex dataset
filename = "data/BrainCellData_268_CellTypes.pkl"

#open the celltypist model
model = celltypist

model = models.Model.load(filename)
