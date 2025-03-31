<p align="center">
	<img src="images/logo.png" width=150>
</p>

<p align="center">
	<a href="https://github.com/ArpiarSaundersLab/tap/tree/main/tests"><img src="https://img.shields.io/badge/build-passing-brightgreen"></a>
	<a href="https://img.shields.io/github/v/release/ArpiarSaundersLab/tap"><img src="https://img.shields.io/github/v/release/ArpiarSaundersLab/tap"></a>
</p>

<br>

# Tropism Analysis Package
**TAP** helps drive discovery of host factors associated with viral tropism in scRNA-seq (or nuclei) datasets.

## Features
- Explore your data in an interactive user interface.
- Uses an ensemble ML approach to rank important features of tropism.
- Easily share results as an HTML file.
- Simple and configurable parameters.

## Usage
TAP accepts an AnnData object and generates an interactive UI in the form of a html file that can be opened in any browser, hosted on a server, or shared. A basic usage example is displayed below. Please see our documentation to go down the rabbit hole of parameters available for tuning your taps.
```python
import tap

parameters = {
	"filename" : "my_file.h5ad",
	"name": "Macaque No Dialout",
	"genes": ["AAV1","AAV2","AAV9","AAV9.Retro"],
	"categories": ["cell_type"],
	"outputPath": ".",
}

results = t.build(parameters)
```

## Live Example

<br>

## Citation

<br>