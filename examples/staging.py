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
