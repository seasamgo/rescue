#' Calculate gini coefficient
#'
#' @param x vector of expression
#'
#' @return gini coefficient
#'

giniCoef <- function(x, weights = rep(1,length(x))) {

  # adapted from R package GiniWegNeg
  # taken from R package Giotto

  dataset = cbind(x, weights)
  ord_x = order(x)
  dataset_ord = dataset[ord_x,]
  x       = dataset_ord[,1]
  weights = dataset_ord[,2]
  N  = sum(weights)
  xw = x*weights
  C_i = cumsum(weights)
  num_1 = sum(xw*C_i)
  num_2 = sum(xw)
  num_3 = sum(xw*weights)
  G_num = (2/N^2)*num_1-(1/N)*num_2-(1/N^2)*num_3
  t_neg = subset(xw, xw<=0)
  T_neg = sum(t_neg)
  T_pos = sum(xw)+abs(T_neg)
  n_RSV = (2*(T_pos+(abs(T_neg)))/N)
  mean_RSV = (n_RSV/2)
  G_RSV = (1/mean_RSV)*G_num
  return(G_RSV)
}


#' Compute Highly Variable Genes
#'
#' @param expression_matrix sparse expression matrix
#' @param method method to calculate highly variable genes
#' @param reverse_log_scale reverse log-scale of expression values
#' @param logbase if reverse_log_scale is TRUE, which log base was used?
#' @param expression_threshold expression threshold to consider a gene detected
#' @param nr_expression_groups number of expression groups for cov_groups
#' @param zscore_threshold zscore to select hvg for cov_groups
#' @param difference_in_variance minimum difference in variance required
#' @return character vector of highly variable genes
#'
#' @export
#'
#' @examples
#' \dontrun{
#' computeHVG(expression_matrix)
#'

computeHVG <- function(
  expression_matrix,
  method = c('cov_groups', 'cov_loess', 'gini_loess'),
  reverse_log_scale = T,
  logbase = exp(1),
  expression_threshold = 0,
  nr_expression_groups = 20,
  zscore_threshold = 1.5,
  difference_in_variance = 1
  ) {

  # method to use
  method = match.arg(method, choices = c('cov_groups', 'cov_loess', 'gini_loess'))

  if(reverse_log_scale == TRUE) {
    expression_matrix = (logbase^expression_matrix)-1
  }

  ## create data.table with relevant statistics ##
  gene_in_cells_detected <- data.table::data.table(
    genes = rownames(expression_matrix),
    nr_cells = Matrix::rowSums(expression_matrix > expression_threshold),
    total_expr = Matrix::rowSums(expression_matrix),
    mean_expr = Matrix::rowMeans(expression_matrix),
    sd = unlist(apply(expression_matrix, 1, sd))
    )
  gene_in_cells_detected[, cov := (sd/mean_expr)]
  gini_level <- unlist(apply(expression_matrix, MARGIN = 1, FUN = function(x) {

    gini = giniCoef(x)
    return(gini)

  }))
  gene_in_cells_detected[, gini := gini_level]

  if(method == 'cov_groups') {

    steps = 1/nr_expression_groups
    prob_sequence = seq(0, 1, steps)
    prob_sequence[length(prob_sequence)] = 1
    expr_group_breaks = stats::quantile(gene_in_cells_detected$mean_expr, probs = prob_sequence)
    expr_groups = cut(x = gene_in_cells_detected$mean_expr, breaks = expr_group_breaks,
                      labels = paste0('group_', 1:nr_expression_groups), include.lowest = T)
    gene_in_cells_detected[, expr_groups := expr_groups]
    gene_in_cells_detected[, cov_group_zscore := scale(cov), by =  expr_groups]
    gene_in_cells_detected[, selected := ifelse(cov_group_zscore > zscore_threshold, 'yes', 'no')]

  } else {

    if(method == 'cov_loess') {
      loess_formula = paste0('cov~mean_expr')
      var_col = 'cov'
    } else if(method == 'gini_loess') {
      loess_formula = paste0('gini~mean_expr')
      var_col = 'gini'
    }

    loess_model_sample <- stats::loess(loess_formula, data = gene_in_cells_detected)
    gene_in_cells_detected$pred_var_genes <- predict(loess_model_sample, newdata = gene_in_cells_detected)
    gene_in_cells_detected[, var_diff := get(var_col)-pred_var_genes, by = 1:nrow(gene_in_cells_detected)]
    setorder(gene_in_cells_detected, -var_diff)
    gene_in_cells_detected[, selected := ifelse(var_diff > difference_in_variance, 'yes', 'no')]
  }

  return(gene_in_cells_detected)

}
