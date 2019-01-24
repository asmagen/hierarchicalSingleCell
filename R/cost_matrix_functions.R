#' Create cost matrix based on distance
#' @param T1_stats A matrix with each row as features for each node
#' @param T2 A matrix with each row as features for each node
#' @param lambda The cost of aligning a node with a empty node
#' @param method The method used to compute correlation (spearman, pearson or euclidean)
#' @return A pairwise distance matrix
#' @export
cost_matrix <- function(T1_stats, T2_stats, lambda = NA, method = 'spearman') {
  cost_matrix <- matrix(NA, nrow = nrow(T1_stats) + 1, ncol = nrow(T2_stats) + 1)
  rownames(cost_matrix) <- c('lambda', rownames(T1_stats))
  colnames(cost_matrix) <- c('lambda', rownames(T2_stats))
  intersect_cols <- intersect(colnames(T1_stats), colnames(T2_stats))
  
  distance_func <- switch(method,
                          spearman  = function(x, y) 1 - cor(x, y, method = 'spearman'),
                          pearson   = function(x, y) 1 - cor(x, y, method = 'pearson'),
                          euclidean = function(x, y) sqrt(sum((x - y) ^ 2)))
  
  for (i in rownames(cost_matrix)) {
    for (j in colnames(cost_matrix)) {
      if (i == 'lambda' & j == 'lambda') cost_matrix[i, j] <- 0
      else if (i == 'lambda' | j == 'lambda') cost_matrix[i, j] <- 1
      else cost_matrix[i, j] <- distance_func(T1_stats[i,intersect_cols], T2_stats[j,intersect_cols])
    }
  }
  if(is.na(lambda)) cost_matrix[is.na(cost_matrix)] <- mean(cost_matrix, na.rm = T)
  return(cost_matrix)
}
