#' create a binary tree node
#' 
#' @param label The label of the node
#' @param left The left tree node
#' @param right The right tree node
#' @return A list representing the tree
#' @export
binary_tree_node <- function(label, left = NULL, right = NULL) {
  node <- list(label = label, left = left, right = right)
  return(node)
}

#' get the children labels
#' 
#' @param node A node object
#' @return A character vector of children's names
get_children <- function(node) {
  children <- c('lambda', 'lambda')
  if (!is.null(node$left)) children[1] <- node$left$label 
  if (!is.null(node$right)) children[2] <- node$right$label 
  return(children)
}

#' Get labels by postorder traversal
#' 
#' @param x A root node of tree
#' @return labels of tree in postorder
#' @export
postorder <- function(x) {
  label_ordered <- c()
  if (!is.null(x$left)) label_ordered <- c(label_ordered, postorder(x$left))
  if (!is.null(x$right)) label_ordered <- c(label_ordered, postorder(x$right))
  label_ordered <- c(label_ordered, x$label)
  return(label_ordered)
}

#' Initialize the matrix to store alignment score
#' @param T1 tree 1
#' @param T2 tree 2
#' @param cost_matrix The cost matrix
#' @return A list of T matrix and F map
#' @export
initialize <- function(T1, T2, cost_matrix) {

  initialize_T1 <- function(x) {
    if (!is.null(x$left)) {
      initialize_T1(x$left)
    }
    if (!is.null(x$right)) {
      initialize_T1(x$right)
    }
    children <- get_children(x)
    F_map[[paste0('1_', x$label)]] <<- sum(T_matrix[children, "lambda"])
    T_matrix[x$label, "lambda"] <<- F_map[[paste0('1_', x$label)]] + cost_matrix[x$label, "lambda"]
  }

  initialize_T2 <- function(x) {
    if (!is.null(x$left)) {
      initialize_T2(x$left)
    }
    if (!is.null(x$right)) {
      initialize_T2(x$right)
    }
    children <- get_children(x)
    F_map[[paste0('2_', x$label)]] <<- sum(T_matrix["lambda", children])
    T_matrix["lambda", x$label] <<- F_map[[paste0('2_', x$label)]] + cost_matrix["lambda", x$label]
  }
  
  ordered_T1 <- postorder(T1)
  ordered_T2 <- postorder(T2)
  T_matrix <- matrix(NA, nrow = length(ordered_T1) + 1, ncol = length(ordered_T2) + 1)
  rownames(T_matrix) <- c('lambda', ordered_T1)
  colnames(T_matrix) <- c('lambda', ordered_T2)
  T_matrix['lambda', 'lambda'] <- 0
  F_map <- hashmap::hashmap('lambda_lambda', 0)
  initialize_T1(T1)
  initialize_T2(T2)
  return(list(T_matrix = T_matrix, F_map = F_map))
}
