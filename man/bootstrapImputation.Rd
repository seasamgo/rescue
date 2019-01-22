\name{bootstrapImputation}
\alias{bootstrapImputation}
\title{Bootstrap Imputation}
\description{
Subsample informative genes, cluster cells using SNN, estimate missing expression values with the distribution mean of means extrapolated from these cell clusterings
}
\usage{
bootstrapImputation(expression_matrix, select_cells = NULL, select_genes = NULL,
  proportion_genes = 0.6, bootstrap_samples = 100, number_pcs = 8, snn_resolution = 0.9,
  impute_index = NULL, use_mclapply = FALSE, cores = 2, return_individual_results = FALSE,
  verbose = FALSE)
}

\arguments{
  \item{expression_matrix}{
A gene by cell log-normalized expression matrix.
}
  \item{select_cells}{
Subset cells if desired.
}
  \item{select_genes}{
A vector of highly variable of differentially expressed gene names, defaults to all genes.
}
  \item{proportion_genes}{
Proportion of informative genes to sample.
}
  \item{bootstrap_samples}{
Number of dimensions to inform SNN clustering.
}
  \item{number_pcs}{
umber of dimensions to inform SNN clustering.
}
  \item{snn_resolution}{
Resolution parameter for SNN.
}
  \item{impute_index}{
Index to impute, will default to all zeroes.
}
  \item{use_mclapply}{
Run in parallel.
}
  \item{cores}{
Number of cores for parallelization.
}
  \item{return_individual_results}{
Return a list of subsampled means.
}
  \item{verbose}{
Print progress output to the console.
}
}
\value{
A list containing the imputed and original expression matrices
}
\examples{
  library(Seurat)
  pbmc_small
  pbmc_small_imputed <- bootstrapImputation(
    expression_matrix = pbmc_small@data,
    select_genes = pbmc_small@var.genes
    )
}

