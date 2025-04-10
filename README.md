<p align="center">
	<img src="assets/images/logo.png" width=150>
</p>

<p align="center">
	<a href="https://github.com/ArpiarSaundersLab/tap/tree/main/tests"><img src="https://img.shields.io/badge/build-passing-brightgreen"></a>
	<a href="https://img.shields.io/github/v/release/ArpiarSaundersLab/tap"><img src="https://img.shields.io/github/v/release/ArpiarSaundersLab/tap"></a>
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

## Usage
TAP accepts an AnnData object or file and generates an interactive UI in the form of a html file that can be opened in any browser, hosted on a server, or shared. A basic usage example is displayed below. Please see our documentation to go down the rabbit hole of parameters available for tuning your taps.
```python
import tap as t

parameters = {
	"filename" : "../assets/data/PCCM.h5ad",
	"name": "PCCM",
	"genes": ["AAV1","AAV2","AAV9"],
	"categories": ["cell_type"],
	"outputPath": ".",
	"outputName": "test.html",
	"runCellTypist" : False,
	"excludeMarkers" : True,
	"mapOnly" : True,
}

results = t.TAP(**parameters)
```

## Live Example

<br>

## Citation

<br>