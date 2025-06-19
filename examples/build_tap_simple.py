import tap as t

parameters = {
	"filename": "filename.h5ad", #Any h5ad file 
	"name": "My first TAP", #This will be displayed in the UI
	"genes": ["AAV1","AAV2","Retro","Anc80"], #viral counts existing as gene/columns in the h5ad file
	"categories": ["cell_type"], #any categorical obs layer field, 
	"categoryNames": ["Cell Type"], #the display name for the category
	"outputPath": ".", #where the output will be saved
	"outputName": "tap.html", #the name of the output file
	"excludeMarkers" : True, #whether to exclude markers of each cluster from the analysis
	"removeOutliers": True, #whether to viral counts that are outliers
	"clusterMethod": "threshold", #the clustering method to use, can be "threshold" or "gaussian"
	"clusterThreshold": 1, #the threshold for clustering, only used if clusterMethod is "threshold"
	"balance": "smote", #if the data is imbalanced, it can be balanced using "smote" or None
}

t.TAP(**parameters)