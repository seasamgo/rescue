
#' Compute Highly Variable Genes
#'
#' @param expression_matrix Expression matrix
#' @param reverse_log_scale Reverse log-scale of expression values
#' @param log_base If reverse_log_scale is TRUE, which log base was used?
#' @param expression_threshold Expression threshold to consider a gene detected
#' @param nr_expression_groups Number of expression groups for cov_groups
#' @param zscore_threshold Z-score to select hvg for cov_groups
#'
#' @return Character vector of highly variable genes
#'
#' @importFrom data.table :=
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
#' ## calculate HVGs
#' hvgs <- computeHVG(m)
#'

computeHVG <- function(
  expression_matrix,
  reverse_log_scale = T,
  log_base = exp(1),
  expression_threshold = 0,
  nr_expression_groups = 20,
  zscore_threshold = 1.5
  ) {

  # data.table variables
  covariance = stdev = mean_expr = expr_groups = cov_group_zscore = selected = NULL

  if(reverse_log_scale == TRUE) {
    expression_matrix = (log_base^expression_matrix)-1
  }

  ## create data.table with relevant statistics ##
  gene_in_cells_detected <- data.table::data.table(
    genes = rownames(expression_matrix),
    nr_cells = Matrix::rowSums(expression_matrix > expression_threshold),
    total_expr = Matrix::rowSums(expression_matrix),
    mean_expr = Matrix::rowMeans(expression_matrix),
    stdev = unlist(apply(expression_matrix, 1, stats::sd))
    )

  gene_in_cells_detected[, covariance := (stdev/mean_expr)]

  ##  bin and score  ##
  steps = 1/nr_expression_groups
  prob_sequence = seq(0, 1, steps)
  prob_sequence[length(prob_sequence)] = 1
  expr_group_breaks = stats::quantile(gene_in_cells_detected$mean_expr, probs = prob_sequence)

  gene_in_cells_detected[, expr_groups := cut(x = gene_in_cells_detected$mean_expr, breaks = expr_group_breaks, labels = paste0('group_', 1:nr_expression_groups), include.lowest = T)]
  gene_in_cells_detected[, cov_group_zscore := scale(covariance), by = expr_groups]
  gene_in_cells_detected[, selected := ifelse(cov_group_zscore > zscore_threshold, 'yes', 'no')]

  return(gene_in_cells_detected)

}
