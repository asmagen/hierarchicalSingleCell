euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

#' Create cost matrix with 1 - cor
#' @param T1_stats A matrix with each row as features for each node
#' @param T2 A matrix with each row as features for each node
#' @param lambda The cost of aligning a node with a empty node
#' @param method The method used to compute correlation (spearman or pearson)
#' @return A pairwise distance matrix (1 - cor)
#' @export
cor_cost_matrix <- function(T1_stats, T2_stats, lambda = NA, method = 'spearman') {
  cost_matrix <- matrix(NA, nrow = nrow(T1_stats) + 1, ncol = nrow(T2_stats) + 1)
  rownames(cost_matrix) <- c('lambda', rownames(T1_stats))
  colnames(cost_matrix) <- c('lambda', rownames(T2_stats))
  intersect_cols <- intersect(colnames(T1_stats), colnames(T2_stats))
  
  for (i in rownames(cost_matrix)) {
    for (j in colnames(cost_matrix)) {
      if (i == 'lambda' & j == 'lambda') cost_matrix[i, j] <- 0
      else if (i == 'lambda' | j == 'lambda') cost_matrix[i, j] <- 1
      else cost_matrix[i, j] <- 1 - cor(T1_stats[i,intersect_cols], T2_stats[j,intersect_cols], method = method)
    }
  }
  if(is.na(lambda)) cost_matrix[is.na(cost_matrix)] <- mean(cost_matrix, na.rm = T)
  return(cost_matrix)
}

#' Create cost matrix with 1 - cor
#' @param T1_stats A matrix with each row as features for each node
#' @param T2 A matrix with each row as features for each node
#' @param lambda The cost of aligning a node with a empty node
#' @return A pairwise distance matrix (1 - cor)
#' @export
euc_cost_matrix <- function(T1_stats, T2_stats, lambda = NA) {
  cost_matrix <- matrix(NA, nrow = nrow(T1_stats) + 1, ncol = nrow(T2_stats) + 1)
  rownames(cost_matrix) <- c('lambda', rownames(T1_stats))
  colnames(cost_matrix) <- c('lambda', rownames(T2_stats))
  intersect_cols <- intersect(colnames(T1_stats), colnames(T2_stats))
  
  for (i in rownames(cost_matrix)) {
    for (j in colnames(cost_matrix)) {
      if (i == 'lambda' & j == 'lambda') cost_matrix[i, j] <- 0
      else if (i == 'lambda' | j == 'lambda') cost_matrix[i, j] <- lambda
      else cost_matrix[i, j] <- euc_dist(T1_stats[i,intersect_cols], T2_stats[j,intersect_cols])
    }
  }
  if(is.na(lambda)) cost_matrix[is.na(cost_matrix)] <- mean(cost_matrix, na.rm = T)
  return(cost_matrix)
}
