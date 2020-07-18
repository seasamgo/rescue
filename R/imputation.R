
#' Sample-mean Estimation
#'
#' Cluster cells using SNN and a list of given genes, estimate missing expression
#' values for each cell-gene combination with the within-cluster non-zero expression
#' mean
#'
#'
#' @param expression_matrix Row by column log-normalized expression matrix
#' @param subset_genes A vector of informative gene names, defaults to all genes
#' @param scale_data Whether to standardize expression by gene, default TRUE
#' @param number_pcs Number of dimensions to inform SNN clustering
#' @param k_neighbors Number of k neighbors to use for NN network
#' @param snn_resolution Resolution parameter for SNN
#' @param impute_index Index to impute, will default to all zeroes
#' @param pseudo_zero Pseudo-zero expression value
#' @param python_path path to your python binary (default = system path)
#' @param verbose Print progress output to the console
#'
#' @return Returns a sparse matrix of class 'dgCMatrix'
#'
#'
#' @export
#'
#' @examples
#' set.seed(0)
#' requireNamespace("Matrix")
#'
#' ## generate (meaningless) counts
#' c1 <- stats::rpois(5e3, 1)
#' c2 <- stats::rpois(5e3, 2)
#' m <- t(
#'   rbind(
#'     matrix(c1, nrow = 20),
#'     matrix(c2, nrow = 20)
#'   )
#' )
#'
#' ## construct an expression matrix m
#' colnames(m) <- paste0('cell', 1:ncol(m))
#' rownames(m) <- paste0('gene', 1:nrow(m))
#' m <- log(m/colSums(m)*1e4 + 1)
#' m <- methods::as(m, 'dgCMatrix')
#'
#' ## impute
#' \donttest{
#' m_imputed <- rescue::sampleImputation(
#'   expression_matrix = m,
#'   k_neighbors = 10
#' )
#' }
#'

sampleImputation <- function(
  expression_matrix,
  subset_genes = NULL,
  scale_data = TRUE,
  number_pcs = 8,
  k_neighbors = 30,
  snn_resolution = .9,
  impute_index = NULL,
  pseudo_zero = NULL,
  python_path = NULL,
  verbose = FALSE
  ) {


  # consider only a subset of genes for clustering:
  # e.g. highly variable genes or random subsets of genes
  if(!is.null(subset_genes)) {
    if(verbose) cat('subsetting genes \n \n')
    my_small_matrix <- expression_matrix[subset_genes, ]
  } else {
    my_small_matrix <- expression_matrix
  }

  # impute a temporary pseudo_zero to ensure SNN convergence
  if(!is.null(pseudo_zero)) my_small_matrix[my_small_matrix == 0] <- pseudo_zero


  ##  STANDARDIZE  ##
  if(verbose) cat('standardizing data \n \n')
  my_small_matrix_scaled <- Matrix::t(scale(Matrix::t(my_small_matrix), center = T, scale = T))


  ##  PCA  ##
  if(verbose) cat('computing principal components \n \n')

  pcs_compute <- min(number_pcs, nrow(x = my_small_matrix_scaled)-1)
  pca_results <- irlba::irlba(A = my_small_matrix_scaled, nv = pcs_compute)
  cell_embeddings <- methods::as(pca_results$v, 'dgCMatrix')
  colnames(cell_embeddings) <- paste0('PC',1:number_pcs)
  rownames(cell_embeddings) <- colnames(expression_matrix)


  ## SNN ##
  if(verbose) cat('constructing nearest neighbors graph \n \n')

  nn <- constructNN(reduced_object = cell_embeddings, k_neighbors = k_neighbors)

  if(verbose) cat('clustering nearest neighbors \n \n')

  cluster_results <- clusterLouvain(nn_network = nn, resolution = snn_resolution, python_path = python_path)

  clusters <- cluster_results$cluster
  names(clusters) <- cluster_results$cell_ID

  ## IMPUTATION step ##
  impute_result_list <- list()

  for(cluster in unique(clusters)) {

    if(verbose) cat('estimating for cluster: ', cluster, '\n')


    # cells that belong to the cluster
    temp_mat <- expression_matrix[, which(colnames(expression_matrix) %in%  names(clusters[clusters == cluster]))]
    temp_nonzero <- temp_mat != 0

    # calculate mean expression and set NaN to 0
    imputed_expr <- Matrix::rowMeans(temp_mat)
    imputed_expr[is.nan(imputed_expr)] <- 0
    imputed_matrix <- Matrix::Matrix(
      data = rep(imputed_expr, ncol(temp_mat)),
      ncol = ncol(temp_mat),
      byrow = FALSE,
      dimnames = list(rownames(temp_mat), colnames(temp_mat))
    )


    # replace zero values (or low values) with imputed values
    if( is.null(impute_index) )
      zero_index <- temp_mat == 0
    else
      zero_index <- impute_index[, which(colnames(expression_matrix) %in%  names(clusters[clusters == cluster]))]


    imp_temp_mat <- temp_mat
    imp_temp_mat[zero_index] <- imputed_matrix[zero_index]
    impute_result_list[[as.character(cluster)]] <- imp_temp_mat


  }


  final_impute_result <- do.call('cbind', impute_result_list)
  final_impute_result <- final_impute_result[match(rownames(expression_matrix), rownames(final_impute_result)),
                                             match(colnames(expression_matrix), colnames(final_impute_result))]


  return(final_impute_result)

}


#' Bootstrap Imputation
#'
#' Subsample informative genes, cluster cells using SNN, estimate missing expression
#' values with the distribution mean of means extrapolated from these cell clusterings
#'
#'
#' @param expression_matrix Row by column log-normalized expression matrix
#' @param select_cells Subset cells if desired
#' @param select_genes A vector of highly variable of differentially expressed gene names,
#' defaults to the most variable
#' @param log_transformed Whether the expression matrix has been log-transformed
#' @param log_base If log-transformed, log-base used
#' @param proportion_genes Proportion of informative genes to sample
#' @param bootstrap_samples Number of samples for the bootstrap
#' @param number_pcs Number of dimensions to inform SNN clustering
#' @param k_neighbors Number of k neighbors to use for NN network
#' @param snn_resolution Resolution parameter for SNN
#' @param impute_index Index to impute, will default to all zeroes
#' @param use_mclapply Run in parallel, default FALSE
#' @param cores Number of cores for parallelization
#' @param return_individual_results Return a list of subsampled means
#' @param python_path path to your python binary (default = system path)
#' @param verbose Print progress output to the console
#'
#' @return Returns a list with the imputed and original expression matrices
#'
#' @export
#'
#' @examples
#' set.seed(0)
#' requireNamespace("Matrix")
#'
#' ## generate (meaningless) counts
#' c1 <- stats::rpois(5e3, 1)
#' c2 <- stats::rpois(5e3, 2)
#' m <- t(
#'   rbind(
#'     matrix(c1, nrow = 20),
#'     matrix(c2, nrow = 20)
#'   )
#' )
#'
#' ## construct an expression matrix m
#' colnames(m) <- paste0('cell', 1:ncol(m))
#' rownames(m) <- paste0('gene', 1:nrow(m))
#' m <- log(m/colSums(m)*1e4 + 1)
#' m <- methods::as(m, 'dgCMatrix')
#'
#' ## impute
#' \donttest{
#' m_imputed <- rescue::bootstrapImputation(
#'   expression_matrix = m,
#'   proportion_genes = .9,
#'   bootstrap_samples = 2,
#'   k_neighbors = 10
#' )
#' }
#'

bootstrapImputation <- function(
  expression_matrix,
  select_cells = NULL,
  select_genes = NULL,
  log_transformed = TRUE,
  log_base = exp(1),
  proportion_genes = 0.6,
  bootstrap_samples = 100,
  number_pcs = 8,
  k_neighbors = 30,
  snn_resolution = .9,
  impute_index = NULL,
  use_mclapply = FALSE,
  cores = 2,
  return_individual_results = FALSE,
  python_path = NULL,
  verbose = FALSE
  ) {

  if(class(expression_matrix) != 'dgCMatrix') expression_matrix <- methods::as(expression_matrix, 'dgCMatrix')

  # test cell subsets
  if(is.null(select_cells)) select_cells = 1:ncol(expression_matrix)
  else if(verbose) cat('subsetting cells \n \n')
  expression_matrix <- expression_matrix[, select_cells]


  ## compute pseudo_zero ##
  pseudo_zero <- min(expression_matrix[which(methods::as(expression_matrix > 0, 'matrix'))])/2


  ## store zero index ##
  if(is.null(impute_index)) impute_index <- expression_matrix == 0


  ## determine variable genes ##
  if(is.null(select_genes)){

    if(verbose) cat('finding variable genes \n \n')

    # data.table variables
    selected = NULL

    hvgs <- computeHVG(expression_matrix, reverse_log_scale = log_transformed, log_base = log_base)
    select_genes <- hvgs[ selected == 'yes', ]$genes

    if(is.null(select_genes)) stop('No HVGs detected by default!')
    else if(verbose) cat('using', length(select_genes), 'variable genes \n \n')
  }


  ## determine number of genes to sample ##
  total_number_of_genes <- length(select_genes)
  number_of_genes_to_use <- floor(total_number_of_genes*proportion_genes)


  ##  BOOTSTRAP  ##
  if(use_mclapply) {
    result_list <- parallel::mclapply(
      X = 1:bootstrap_samples,
      mc.preschedule = FALSE,
      mc.cores = cores,

      FUN = function(round) {

      gene_sample <- sample(x = 1:total_number_of_genes, size = number_of_genes_to_use, replace = F)
      genes_to_use <- select_genes[gene_sample]

      temp_impute <- sampleImputation(
        expression_matrix = expression_matrix,
        subset_genes = genes_to_use,
        number_pcs = number_pcs,
        k_neighbors = k_neighbors,
        snn_resolution = snn_resolution,
        impute_index = impute_index,
        pseudo_zero = pseudo_zero,
        python_path = python_path,
        verbose = verbose
      )

      return(temp_impute)

    })
  } else {
    result_list <- list()
    for(round in 1:bootstrap_samples) {

      gene_sample  = sample(x = 1:total_number_of_genes, size = number_of_genes_to_use, replace = F)
      genes_to_use = select_genes[gene_sample]

      temp_impute = sampleImputation(
        expression_matrix = expression_matrix,
        subset_genes = genes_to_use,
        number_pcs = number_pcs,
        k_neighbors = k_neighbors,
        snn_resolution = snn_resolution,
        impute_index = impute_index,
        pseudo_zero = pseudo_zero,
        python_path = python_path,
        verbose = verbose
      )

      result_list[[as.character(round)]] <- temp_impute
      if(verbose) cat('sample: ', round, '\n \n')
    }
  }


  ## calculate average of bootstrapped sample means ##
  final_imputation <- Reduce('+', result_list)/bootstrap_samples
  final_imputation <- final_imputation


  ## return data ##
  if(return_individual_results) {
    return(
      list(
        individual_results = result_list,
        final_imputation = final_imputation,
        original_matrix = expression_matrix
      )
    )
  } else {
    return(
      list(
        final_imputation = final_imputation,
        original_matrix = expression_matrix
      )
    )
  }

}


