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
devtools::install_github("wang-lab/SpatialDeX", build_vignettes = TRUE)
```


## Running the example
An example of raw UMI matrix from a subset of a tumor sample sequenced by 10X Genomics platform is included with this package named, example_raw_count.

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

## Running deconvolution
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

***lambda*** lambda is the hyperparameter used for regularization in ridge regression. Default is 1.

***sam.name*** is sample name.

***LOW.DR*** is minimal population fractions of genes for smoothing for CNV detection. Default is 0.03.

***UP.DR*** is minimal population fractions of genes for segmentation for CNV detection. Default is 0.2.

***ngene.chr*** is minimal number of genes per chromosome for cell filtering for CNV detection. Default is 5.

***KS.cut*** is segmentation parameters, input 0 to 1; larger looser criteria for CNV detection. Default is 0.1.

***win.size*** is minimal window sizes for segmentation for CNV detection. Default is 25.

***top_cell_type*** define the number of cell types to retain based on the proportion of each cell type. Choose an integer from 0 to 9.


### Outputs

Result is a dataframe with the following fields:  
1. spot barcode  
2. cell type proportions for individual major cell types

## Get the weighted feature coefficient by users

This function allows users to obtain their own weighted feature coefficient matrix. Below provides an example of identity of cell types across different spots. This identify matrix can be provided by users. Each row represents a spot, and each column represents a cell type. A value of 1 indicates the presence of that cell type in the spot, while 0 indicates its absence.


```
data(example_raw_count)
raw_exprs <- example_raw_count[,1:10]
cell_identity <- matrix(0, nrow = 10, ncol = 9)
rownames(cell_identity) <- colnames(raw_exprs)
colnames(cell_identity) <- paste0("cell_type", 1:9)
for (i in 1:10) {
  n_ones <- sample(0:3, 1)
  if (n_ones > 0) {
    ones_positions <- sample(1:9, n_ones)
    cell_identity[i, ones_positions] <- 1
  }
}
```

Run the code below to get the weighted feature coefficient matrix.

```
feat_coef_mat<- get_feat_coefficient(cell_identity,
                                raw_exprs=raw_exprs,
                                LAMBDA=0.01,
                                n_permutation=10000,
                                pvalue=0.05,
                                sam.name="sam.name",
                                LOW.DR=0.03,
                                UP.DR=0.2,
                                ngene.chr=5,
                                KS.cut=0.1,
                                win.size=25,
                                seed=1234)
```

### Inputs

***raw_exprs*** is raw gene expression matrix, with gene symbols in rows and spot names in columns.

***cell_identity*** is cell type identity matrix. Each row represents a spot, and each column represents a cell type. A value of 1 indicates the presence of that cell type in the spot, while 0 indicates its absence.

### Parameters

***n_permutation*** n_permutation An integer indicating the number of permutations for randomly shuffling the cell labels of the input matrix to create a null distribution. Default is 10,000.

***LAMBDA*** A numeric value for the regularization parameter. Default is 0.01.

***sam.name*** Sample name

***LOW.DR*** Minimal population fractions of genes for smoothing for CNV detection. Default is 0.03.

***UP.DR*** Minimal population fractions of genes for segmentation for CNV detection. Default is 0.2.

***ngene.chr*** Minimal number of genes per chromosome for cell filtering for CNV detection. Default is 5.

***KS.cut*** Segmentation parameters, input 0 to 1; larger looser criteria for CNV detection. Default is 0.1.

***win.size*** Minimal window sizes for segmentation for CNV detection. Default is 25.

***pvalue*** Significance threshold for weighted feature coefficients. Calculated by permutation.

***seed*** Seed for random number generation for reproducibility.

### Outputs

Result is a matrix of coefficients for each feature to each cell type is generated. Each cell in the matrix contains a coefficient indicating the strength of the relationship between a feature and a given cell type.
