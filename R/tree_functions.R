#' create a binary tree node
#' 
#' @param label The label of the node
#' @param left The left tree node
#' @param right The right tree node
#' @return A list representing the tree
#' @export
binary_tree <- function(label, left = NULL, right = NULL) {
  tree <- list(label = label, left = left, right = right)
  class(tree) <- 'binary_tree'
  return(tree)
}

make_children_map <- function()

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

#' A skeleton function
postorder <- function(x, f) {
  if (!is.null(x$left)) postorder(x$left, f)
  if (!is.null(x$right)) postorder(x$right, f)
  f(x)
}

#' Get labels by postorder traversal
#' 
#' @param x A root node of tree
#' @return labels of tree in postorder
#' @export
postorder_labels <- function(x) {
  labels <- c()
  postorder(x, function(t) labels <<- c(labels, t$label))
  return(labels)
}

#' Create alignment object
#' @param T1 tree 1
#' @param T2 tree 2
#' @param cost_matrix The cost matrix
#' @return An alignment object
#' @export
create_align_object <- function(T1, T2, cost_matrix) {
  ordered_T1 <- postorder_labels(T1)
  ordered_T2 <- postorder_labels(T2)
  T_matrix <- matrix(NA, nrow = length(ordered_T1) + 1, ncol = length(ordered_T2) + 1)
  rownames(T_matrix) <- c('lambda', ordered_T1)
  colnames(T_matrix) <- c('lambda', ordered_T2)
  T_matrix['lambda', 'lambda'] <- 0
  F_map <- hashmap::hashmap('lambda_lambda', 0)
  list(T1 = T1, T2 = T2, cost_matrix = cost_matrix, T_matrix = T_matrix, F_map = F_map)
}

initialize_T1 <- function(align_object) {
  postorder(align_object$T1,
            function(t) {
              children <- get_children(t)
              align_object$F_map[[paste_collapse(1, t$label, 1, 2, 'lambda')]] <<- sum(align_object$T_matrix[children, "lambda"])
              align_object$T_matrix[t$label, "lambda"] <<- align_object$F_map[[paste_collapse(1, t$label, 1, 2, 'lambda')]] + align_object$cost_matrix[t$label, "lambda"]
            })
  return(align_object)
}

initialize_T2 <- function(align_object) {
  postorder(align_object$T2,
    function(t) {
      children <- get_children(t)
      align_object$F_map[[paste_collapse('lambda', 2, t$label, 1, 2)]] <<- sum(align_object$T_matrix["lambda", children])
      align_object$T_matrix["lambda", t$label] <<- align_object$F_map[[paste_collapse('lambda', 2, t$label, 1, 2)]] + align_object$cost_matrix["lambda", t$label]
    }
  )
  return(align_object)
}

initialize <- function(align_object) {
  align_object <- initialize_T1(align_object)
  align_object <- initialize_T2(align_object)
  align_object$T1_children <- create_children_map(align_object$T1)
  align_object$T2_children <- create_children_map(align_object$T2)
  return(align_object)
}

create_children_map <- function(tree) {
  childrens <- c()
  labels <- c()
  postorder(tree, function(x) {
    left_children <- ifelse(is.null(x$left), 'lambda', x$left$label)
    right_children <- ifelse(is.null(x$right), 'lambda', x$right$label)
    childrens <<- rbind(childrens, c(left_children, right_children))
    labels <<- c(labels, x$label)
  })
  rownames(childrens) <- labels
  colnames(childrens) <- c('left', 'right')
  return(childrens)
}

lambda_T2 <- function(align_obj, i, j) {
  subtree_cost <- c()
  T2_children_j <- align_obj$T2_children[j, ]
  for (r in T2_children_j) {
    subtree_cost <- c(subtree_cost, align_obj$T_matrix[i, r] - align_obj$T_matrix['lambda', r])
  }
  align_obj$T_matrix['lambda', j] + min(subtree_cost)
}

T1_lambda <- function(align_obj, i, j) {
  subtree_cost <- c()
  T1_children_i <- align_obj$T1_children[i, ]
  for (r in T1_children_i) {
    subtree_cost <- c(subtree_cost, align_obj$T_matrix[r, j] - align_obj$T_matrix[r, 'lambda'])
  }
  align_obj$T_matrix[i, 'lambda'] + min(subtree_cost)
}

T1_T2 <- function(align_obj, i, j) {
  align_obj$F_map[[paste(1, i, 1, 2, 2, j, 1, 2, collapse = '_')]] + align_obj$cost_matrix[i, j]
}

paste_collapse <- function(...) paste0(c(...), collapse = '_')

fill_matrix <- function(align_obj) {
  align_obj$choice <- align_obj$T_matrix
  align_obj$choice[,] <- NA
  for (i in rownames(align_obj$T_matrix)[-1]) {
    cost <- c()
    for (j in colnames(align_obj$T_matrix)[-1]) {
      
      for (s in seq(2)) {
        align_obj <- fill_F(align_obj, paste_collapse(1, i, s, 2), 
                            paste_collapse(2, j, 1, 2))
      }
      
      for (t in seq(2)) {
        align_obj <- fill_F(align_obj, paste_collapse(1, i, 1, 2), 
                            paste_collapse(2, j, t, 2))
      }
      
      cost <- c(cost, lambda_T2(align_obj, i, j))
      cost <- c(cost, T1_lambda(align_obj, i, j))
      cost <- c(cost, T1_T2(align_obj, i, j))
      align_obj$choice[i, j] <- which.min(cost)
      align_obj$T_matrix[i, j] <- min(cost)
    }
  }
  return(align_obj)
}

fill_F <- function(align_obj, x, y) {
  x <- strsplit(x, split = '_')[[1]][-1]
  i <- x[1]
  s <- as.numeric(x[2])
  mi <- as.numeric(x[3])
  y <- strsplit(y, split = '_')[[1]][-1]
  j <- y[1]
  t <- as.numeric(y[2])
  nj <- as.numeric(y[3])
  #browser()
  align_obj$F_map[[paste_collapse(1, i, s, s-1, 2, j, t, t - 1)]] <- 0
  for (p in seq(s, mi)) {
    align_obj$F_map[[paste_collapse(1, i, s, p, 2, j, t, t - 1)]] <- 
      align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, t - 1)]] +
      align_obj$T_matrix[align_obj$T1_children[i, p], 'lambda']
  }
  
  for (q in seq(t, nj)) {
    align_obj$F_map[[paste_collapse(1, i, s, s - 1, 2, j, t, q)]] <- 
      align_obj$F_map[[paste_collapse(1, i, s, s - 1, 2, j, t, q - 1)]] +
      align_obj$T_matrix['lambda', align_obj$T2_children[j, q]]
  }
  
  for (p in seq(s, mi)) {
    for (q in seq(t, nj)) {
      align_obj <- fill_F_helper(align_obj, i, s, p, j, t, q)
    }
  }
  return(align_obj)
}

fill_F_helper <- function(align_obj, i, s, p, j, t, q) {
  case_1 <- align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, q)]] +
    align_obj$T_matrix[align_obj$T1_children[i, p], 'lambda']
  case_2 <- align_obj$F_map[[paste_collapse(1, i, s, p, 2, j, t, q - 1)]] +
    align_obj$T_matrix['lambda', align_obj$T2_children[j, q]]
  case_3 <- align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, q - 1)]] +
    align_obj$T_matrix[align_obj$T1_children[i, p], align_obj$T2_children[j, q]]
  case_4_temp <- c()
  for (k in seq(s, p)) {
    case_4_temp <- c(case_4_temp, align_obj$F_map[[paste_collapse(1, i, s, k - 1, 2, j, t, q - 1)]] +
                       align_obj$F_map[[paste_collapse(1, i, k, p, 2, align_obj$T2_children[j, q], 1, 2)]])
  }
  case_4 <- align_obj$cost_matrix['lambda', align_obj$T2_children[j, q]] + min(case_4_temp)
  
  case_5_temp <- c()
  for (k in seq(t, q)) {
    case_5_temp <- c(case_5_temp, align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, k - 1)]] +
                       align_obj$F_map[[paste_collapse(1, i, align_obj$T1_children[i, p], 1, 2, 2, j, k, q)]])
  }
  case_5 <- align_obj$cost_matrix[align_obj$T1_children[i, p], 'lambda'] + min(case_5_temp)
  
  align_obj$F_map[[paste_collapse(1, i, s, p, 2, j, t, q)]] <- min(case_1, case_2, case_3, case_4, case_5)
  return(align_obj)
}