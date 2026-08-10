<p align="center">
	<img src="assets/images/logo.png" width=150>
</p>

<p align="center">
	<a href="https://github.com/ArpiarSaundersLab/tap/tree/main/tests"><img src="https://img.shields.io/badge/build-passing-brightgreen"></a>
</p>

<br>

# Tropism Analysis Package
**TAP** helps drive discovery of host factors associated with viral tropism in sc/nRNA-seq datasets by automating infectivity clustering, then applying an essemble feature selection pipeline to each subset of cells. The results of each subset are stored and presented in a easy-to-use intuitive interface.

## Features
- Explore your data in an interactive user interface.
- Uses an ensemble ML approach to rank important features of tropism.
- Easily share results as an HTML file.
- Simple and configurable parameters.

## Installation
We recommend installing in a fresh python (>=3.12) environment.
```
pip install scTap
```

## Live Demo
A curated selection of multiplexed AAV TAPs to explore.<br />
<a href="https://www.aavdb.com" target="_blank">aavdb.com</a>

## Example Notebooks
Notebook and starter examples are organized in the examples area of this repository.

- Browse the examples directory: [examples](examples)
- Minimal TAP script: [examples/build_tap_simple.py](examples/build_tap_simple.py)
- Notebook example: [examples/full_example.ipynb](examples/full_example.ipynb)


## Usage
TAP accepts an AnnData object or file and generates an interactive HTML report that can be opened in any browser, hosted on a server, or shared. A basic example is shown below. The package expects raw counts in the AnnData object, so the input data should include the original expression values in `adata.raw.X` for the most reliable results.

```python
import tap as t

parameters = {
    "filename": "../assets/data/PCCM.h5ad",
    "name": "PCCM",
    "genes": ["AAV1", "AAV2", "AAV9"],
    "categories": ["cell_type"],
    "outputPath": ".",
    "outputName": "test.html",
    "runCellTypist": False,
    "excludeMarkers": True,
    "mapOnly": True,
}

results = t.TAP(**parameters)
```

### TAP constructor parameters
The `TAP` class accepts the following constructor arguments. These options control input data handling, clustering, feature selection, plotting, and output generation.

| Parameter | Default | Description |
| --- | --- | --- |
| `name` | `"Heatmap"` | Name used for the report and related metadata. |
| `filename` | `None` | Path to an AnnData file on disk. |
| `adataObject` | `None` | An AnnData object to use directly instead of loading from `filename`. |
| `annSqlDB` | `None` | Connection or identifier for an AnnSQL database to convert into AnnData. |
| `categories` | `[]` | One or more metadata columns used to define the grouping labels for analysis. |
| `categoryNames` | `[]` | Optional display names for the categories used in the report. |
| `genes` | `[]` | Gene names or features of interest to include when building the analysis. |
| `exclude` | `[]` | Genes or features to exclude from the analysis. |
| `useRaw` | `False` | Whether to use raw values from the AnnData object when available. |
| `outputPath` | `"."` | Directory where the output HTML file and related files will be written. |
| `outputName` | `"tap.html"` | Name of the generated HTML file. |
| `cpus` | `2` | Number of CPU cores to use for parallel work where supported. |
| `mapOnly` | `False` | If `True`, generate the mapping outputs without the full reporting workflow. |
| `showDE` | `True` | Whether to include differential expression summaries in the output. |
| `showRF` | `True` | Whether to include random forest results in the output. |
| `showUMAP` | `True` | Whether to include UMAP visualizations in the output. |
| `showDetails` | `True` | Whether to include detailed per-cluster or per-feature information. |
| `clusterMethod` | `"threshold"` | Method used to define the clustering threshold for analysis. |
| `clusterThreshold` | `1` | Threshold value used by the clustering procedure. |
| `useLog10` | `False` | Whether to log transform values before downstream analysis. |
| `showPlots` | `False` | Whether to include additional plots in the output. |
| `excludeGenes` | `[]` | Genes to remove from the analysis before modeling. |
| `deMethod` | `"wilcoxon"` | Statistical method used for differential expression testing. |
| `rfHyperParameterTune` | `False` | Whether to tune random forest hyperparameters. |
| `rfHyperParameterIterations` | `50` | Number of iterations used for random forest hyperparameter tuning. |
| `rfType` | `"classifier"` | Type of random forest model to use, either `classifier` or `regressor`. |
| `balance` | `None` | Optional balancing strategy for classification tasks. |
| `useAllGenes` | `False` | If `True`, use all genes from the dataset instead of only the genes present in the input. |
| `rfPermuteFeatureImportance` | `False` | Whether to estimate feature importance using permutation testing. |
| `rfPermuteRepeats` | `2` | Number of repeats used for permutation-based feature importance. |
| `runCellTypist` | `False` | Whether to run CellTypist predictions as part of the workflow. |
| `cellTypistModel` | `"Mouse_Whole_Brain.pkl"` | CellTypist model file to use when `runCellTypist` is enabled. |
| `cellTypistPlots` | `False` | Whether to generate plots for CellTypist results. |
| `cellTypistLevels` | `3` | Number of CellTypist prediction levels to include. |
| `minify` | `False` | Whether to minify the generated HTML output. |
| `excludeMarkers` | `False` | Whether to exclude marker genes based on the marker determination step. |
| `removeOutliers` | `True` | Whether to remove outlier cells before analysis. |
| `minCells` | `20` | Minimum number of cells required for a group to be retained. |
| `minSeroTypeCells` | `20` | Minimum number of cells required for a serotype group. |
| `clusterThresholdGreaterThanOrEqual` | `None` | Optional lower bound for the clustering threshold comparison. |
| `clusterThresholdLessThanOrEqual` | `None` | Optional upper bound for the clustering threshold comparison. |
| `remove_outlier_upper_percentile` | `99` | Upper percentile cutoff used when removing outlier cells. |
| `remove_outlier_lower_percentile` | `0` | Lower percentile cutoff used when removing outlier cells. |
| `replicate_mode` | `False` | Enables replicate-oriented processing behavior. |

Commonly used options include `genes`, `categories`, `outputPath`, `outputName`, `excludeMarkers`, `mapOnly`, and `runCellTypist`.


## Citation
Pavan Kenny et al. Tropism Analysis Package: Interactive Machine Learning Software to Identify Viral Host Factors Through Single-Cell Host-Virus mRNA Profiling. bioRxiv. Year; DOI: [add DOI].
<br>