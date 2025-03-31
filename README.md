<p align="center">
	<img src="images/logo.png" width=150>
</p>

<p align="center">
	<a href="https://github.com/ArpiarSaundersLab/tap/tree/main/tests"><img src="https://img.shields.io/badge/build-passing-brightgreen"></a>
	<a href="https://img.shields.io/github/v/release/ArpiarSaundersLab/tap"><img src="https://img.shields.io/github/v/release/ArpiarSaundersLab/tap"></a>
</p>

<br>

# Tropism Analysis Package
*TAP* helps drive discovery of host factors associated with viral tropism in scRNA-seq (or nuclei) datasets.

## Features
- Explore your data in an interactive user interface.
- Uses an ensemble ML approach to rank important features of tropism

<br>

## Usage
TAP accepts an AnnData object and generates an interactive UI in the form of a html file that can be opened in any browser locally, hosted on a server, or shared. 
```python
import tap

parameters = {
	"filename" : "my_file.h5ad"
	"genes": ["AAV1","AAV2","AAV9","AAV9.Retro"]
	"categories" ["cell_type"]
}

results = t.build(adata,parameters)
```

<br>

## Live Example

<br>
<br>

## Citation
<br>
<br>