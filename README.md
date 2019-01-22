
<!-- README.md is generated from README.Rmd. Please edit that file -->
README
------

This package provides a bootstrap imputation method for dropout events in single-cell RNA-seq data.

Installation
------------

``` r
install.packages("devtools", repos="http://cran.rstudio.com/")
library(devtools)
devtools::install_github("seasamgo/rescue")
library(rescue)
```

Method
------

`bootstrapImputation` takes a log-normalized expression matrix and returns a list containing the imputed and original matrices.

``` r
bootstrapImputation(
  expression_matrix,                  # expression matrix
  select_cells = NULL,                # subset cells
  select_genes = NULL,                # informative genes
  proportion_genes = 0.6,             # proportion of genes to sample
  bootstrap_samples = 100,            # number of samples
  number_pcs = 8,                     # number of PC's to consider
  snn_resolution = 0.9,               # clustering resolution
  impute_index = NULL,                # specify non-default gene-cell counts
  use_mclapply = FALSE,               # run in parallel
  cores = 2,                          # number of parallel cores
  return_individual_results = FALSE,  # return sample means
  verbose = FALSE                     # print progress to console
  )
```

Names of informative genes (e.g. highly variable or differentially expressed genes) may be indicated with `select_genes`. For more, please view the help files.
