\name{sampleImputation}
\alias{sampleImputation}
\title{
Sample-mean Estimation
}
\description{
Cluster cells using SNN and a list of given genes, estimate missing expression values for each cell-gene combination with the within-cluster non-zero expression mean
}
\usage{
sampleImputation(expression_matrix, subset_genes = NULL, scale_data = TRUE,
  number_pcs = 8, snn_resolution = 0.9, impute_index = NULL, pseudo_zero = NULL,
  verbose = FALSE)
}
\arguments{
  \item{expression_matrix}{
Gene by cell log-normalized expression matrix.
}
  \item{subset_genes}{
A vector of informative gene names, defaults to all genes.
}
  \item{scale_data}{
Whether to standardize expression by gene, default TRUE.
}
  \item{number_pcs}{
Number of dimensions to inform SNN clustering.
}
  \item{snn_resolution}{
Resolution parameter for SNN.
}
  \item{impute_index}{
Index to impute, will default to all zeroes.
}
  \item{pseudo_zero}{
Pseudo-zero expression value.
}
  \item{verbose}{
Print progress output to the console, default FALSE.
}
}
\value{
A sparse matrix of class 'dgCMatrix'.
}
\examples{
library(Seurat)
pbmc_small
pbmc_small_imputed <- sampleImputation(
  expression_matrix = pbmc_small@data,
  subset_genes = pbmc_small@var.genes
  )
}
