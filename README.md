# SpatialDeX
The rapid development of spatial transcriptomics (ST) technologies has enabled transcriptome-wide profiling of gene expression in tissue sections. Despite the emergence of single-cell resolution platforms, most ST sequencing studies still operate at the multi-cell resolution. Consequently, deconvolution of cell identities within the spatial spots has become imperative for characterizing cell type-specific spatial organization.  To this end, we introduce SpatialDeX, **Spatial** **D**econvolution **Ex**plorer, a regression model-based method for estimating cell type proportions in tumor ST spots.  

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for prediction and testing purposes.

### Prerequisites
The following components were used to run the package. Different versions of each may or may not be compatible.

- R (>= 3.5.0)
- GSVA 
- copykat: copykat[V1.1.0] (https://github.com/navinlabcode/copykat)
- dplyr

### Installing
Installing SpatialDeX from GitHub


```
library(devtools)
install_github("wang-lab/SpatialDeX")
```


## Running the example
An example raw UMI matrix from a subset of a tumor sample sequenced by 10X Genomics platform is included with this package named, example_raw_count.

To test the package, run this line of code in R/Rstudio:

```
library(SpatialDeX)
data(example_raw_count)
res<- SpatialDeX(raw_exprs=example_raw_count,
                      lambda=1,
                      sam.name="test",
                      LOW.DR=0.03,
                      UP.DR=0.2,
                      ngene.chr=5,
                      KS.cut=0.1,
                      win.size=25,
                      top_cell_type=7)
```

## Running SpatialDeX
```
res<- SpatialDeX(raw_exprs,
                      lambda,
                      sam.name,
                      LOW.DR,
                      UP.DR,
                      ngene.chr,
                      KS.cut,
                      win.size,
                      top_cell_type=NULL)

```
### Inputs

***raw_exprs*** is raw gene expression matrix, with gene symbols in rows and spot names in columns.

### Parameters
***lambda*** lambda is the hyperparameter used for regularization in ridge regression.

***sam.name*** is sample name.

***LOW.DR*** is minimal population fractions of genes for smoothing for copykat.

***UP.DR*** is minimal population fractions of genes for segmentation for copykat.

***ngene.chr*** is minimal number of genes per chromosome for cell filtering for copykat.

***KS.cut*** is segmentation parameters, input 0 to 1; larger looser criteria for copykat.

***win.size*** is minimal window sizes for segmentation for copykat.

***top_cell_type*** defind the number of cell types to retain based on the proportion of each cell type. Choose an integer from 0 to 9.



### Outputs
Result is a dataframe with the following fields:  
1. spot barcode  
2. cell type proportions for individual major cell types


