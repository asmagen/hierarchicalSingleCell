paste_collapse <- function(...) paste0(c(...), collapse = ' ')

#' Align two binary ordered trees
#' @param T1 A binary tree
#' @param T2 A binary tree
#' @param cost_matrix Cost matrix between each pair of nodes in T1 and T2
#' @return An alignment object
#' @export
align2 <- function(T1, T2, cost_matrix) {
  align_obj <- create_align_object2(T1, T2, cost_matrix)
  align_obj <- fill_matrix2(align_obj)
  align_obj <- traceback(align_obj)
  align_obj <- build_tree(align_obj)
  return(align_obj)
}


#' Create alignment object
#' @param T1 tree 1
#' @param T2 tree 2
#' @param cost_matrix The cost matrix
#' @return An alignment object
#' @export
create_align_object2 <- function(T1, T2, cost_matrix) {
  ordered_T1 <- postorder_labels2(T1)
  ordered_T2 <- postorder_labels2(T2)
  T_matrix <- matrix(NA, nrow = length(ordered_T1) + 1, ncol = length(ordered_T2) + 1)
  rownames(T_matrix) <- c('lambda', ordered_T1)
  colnames(T_matrix) <- c('lambda', ordered_T2)
  T_matrix['lambda', 'lambda'] <- 0
  F_map <- hashmap::hashmap('lambda lambda', 0)
  T1_children <- create_children_map2(T1)
  T2_children <- create_children_map2(T2)
  align_object <- list(T1 = T1, T2 = T2, cost_matrix = cost_matrix, T_matrix = T_matrix, F_map = F_map,
                    T1_children = T1_children, T2_children = T2_children)
  align_object <- initialize_T12(align_object)
  align_object <- initialize_T22(align_object)
  return(align_object)
}

create_children_map2 <- function(tree) {
  childrens <- list()
  postorder2(tree, function(x) {
    if (!is.null(x$children)) {
      childrens[[x$label]] <<- map_chr(x$children, "label")
    } else {
      childrens[[x$label]] <<- "lambda"
    }
  })
  return(childrens)
}

initialize_T12 <- function(align_object) {
  T1_node <- rownames(align_object$T_matrix)[-1]
  for (node in T1_node) {
    children <- align_object$T1_children[[node]]
    F_node_lambda <- align_object$T_matrix[children, "lambda"]
    n_children <- length(children)
    for (i in seq(n_children)) {
      for (j in seq(i, n_children)) {
      align_object$F_map[[paste_collapse(1, node, i, j, 'lambda')]] <- sum(F_node_lambda[i:j])
      }
    }
    align_object$T_matrix[node, "lambda"] <- sum(F_node_lambda) + align_object$cost_matrix[node, "lambda"]
  }
  return(align_object)
}

initialize_T22 <- function(align_object) {
  T2_node <- colnames(align_object$T_matrix)[-1]
  for (node in T2_node) {
    children <- align_object$T2_children[[node]]
    lambda_F_node <- align_object$T_matrix["lambda", children]
    n_children <- length(children)
    for (i in seq(n_children)) {
      for (j in seq(i, n_children)) {
        align_object$F_map[[paste_collapse('lambda', 2, node, i, j)]] <- sum(lambda_F_node[i:j])
      }
    }
    align_object$T_matrix["lambda", node] <- sum(lambda_F_node) + align_object$cost_matrix["lambda", node]
  }
  return(align_object)
}

#' Fill T matrix and F map
#' 
#' @param align_obj An alignment object
#' @return An alignment object with T matrix and F map filled
#' @export
fill_matrix2 <- function(align_obj) {
  align_obj$T_choice <- align_obj$T_matrix
  align_obj$T_choice[,] <- NA
  align_obj$F_choice <- list()
  
  for (i in rownames(align_obj$T_matrix)[-1]) {
    for (j in colnames(align_obj$T_matrix)[-1]) {
      mi <- length(align_obj$T1_children[[i]])
      nj <- length(align_obj$T2_children[[j]])
      cost <- c()
      for (s in seq(mi)) {
        align_obj <- fill_F2(align_obj, paste_collapse(1, i, s, mi), 
                            paste_collapse(2, j, 1, nj))
      }
      
      for (t in seq(nj)) {
        align_obj <- fill_F2(align_obj, paste_collapse(1, i, 1, mi), 
                            paste_collapse(2, j, t, nj))
      }
      
      cost <- c(cost, lambda_T22(align_obj, i, j))
      cost <- c(cost, T1_lambda2(align_obj, i, j))
      cost <- c(cost, T1_T22(align_obj, i, j))
      # print(cost)
      align_obj$T_choice[i, j] <- which.min(cost)
      # print(paste0(c(i, j, cost), collapse = ', '))
      align_obj$T_matrix[i, j] <- min(cost, na.rm = T)
    }
  }
  return(align_obj)
}

fill_F2 <- function(align_obj, x, y) {
  x <- strsplit(x, split = ' ')[[1]][-1]
  i <- x[1]
  s <- as.numeric(x[2])
  mi <- as.numeric(x[3])
  y <- strsplit(y, split = ' ')[[1]][-1]
  j <- y[1]
  t <- as.numeric(y[2])
  nj <- as.numeric(y[3])
  align_obj$F_map[[paste_collapse(1, i, s, s-1, 2, j, t, t - 1)]] <- 0
  for (p in seq(s, mi)) {
    # print(align_obj$T1_children[[i]][p])
    align_obj$F_map[[paste_collapse(1, i, s, p, 2, j, t, t - 1)]] <- 
      align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, t - 1)]] +
      align_obj$T_matrix[align_obj$T1_children[[i]][p], 'lambda']
  }
  
  for (q in seq(t, nj)) {
    align_obj$F_map[[paste_collapse(1, i, s, s - 1, 2, j, t, q)]] <- 
      align_obj$F_map[[paste_collapse(1, i, s, s - 1, 2, j, t, q - 1)]] +
      align_obj$T_matrix['lambda', align_obj$T2_children[[j]][q]]
  }
  
  for (p in seq(s, mi)) {
    for (q in seq(t, nj)) {
      align_obj <- fill_F_helper2(align_obj, i, s, p, j, t, q)
    }
  }
  return(align_obj)
}

fill_F_helper2 <- function(align_obj, i, s, p, j, t, q) {
  case_1 <- align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, q)]] +
    align_obj$T_matrix[align_obj$T1_children[[i]][p], 'lambda']
  case_2 <- align_obj$F_map[[paste_collapse(1, i, s, p, 2, j, t, q - 1)]] +
    align_obj$T_matrix['lambda', align_obj$T2_children[[j]][q]]
  case_3 <- align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, q - 1)]] +
    align_obj$T_matrix[align_obj$T1_children[[i]][p], align_obj$T2_children[[j]][q]]
  case_4 <- fill_case_42(align_obj, i, s, p, j, t, q)
  case_5 <- fill_case_52(align_obj, i, s, p, j, t, q)
  
  all <- c(case_1, case_2, case_3, case_4, case_5)
  loc <- paste_collapse(1, i, s, p, 2, j, t, q)
  # print(paste0(loc, ":", all))
  align_obj$F_map[[loc]] <- min(all)
  align_obj$F_choice[[loc]] <- which.min(all)
  return(align_obj)
}

fill_case_42 <- function(align_obj, i, s, p, j, t, q) {
  case_4 <- c()
  for (k in seq(s, p)) {
    idx1 <- paste_collapse(1, i, k, p)
    T2_j_q <- align_obj$T2_children[[j]][q]
    if (T2_j_q == 'lambda') idx2 <- 'lambda'
    else idx2 <- paste_collapse(2, T2_j_q, 1, length(align_obj$T2_children[[T2_j_q]]))
    case_4 <- c(case_4, align_obj$F_map[[paste_collapse(1, i, s, k - 1, 2, j, t, q - 1)]] + align_obj$F_map[[paste_collapse(idx1, idx2)]])
  }
  case_4 <- align_obj$cost_matrix['lambda', align_obj$T2_children[[j]][q]] + min(case_4)
  return(case_4)
}

fill_case_52 <- function(align_obj, i, s, p, j, t, q) {
  case_5 <- c()
  for (k in seq(t, q)) {
    T1_i_p <- align_obj$T1_children[[i]][p]
    if (T1_i_p == 'lambda') idx1 <- 'lambda'
    else idx1 <- paste_collapse(1, T1_i_p, 1, length(align_obj$T1_children[[T1_i_p]]))
    idx2 <- paste_collapse(2, j, k, q)
    # print(paste0(paste_collapse(idx1, idx2), align_obj$F_map[[paste_collapse(idx1, idx2)]]))
    case_5 <- c(case_5, align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, k - 1)]] + align_obj$F_map[[paste_collapse(idx1, idx2)]])
  }
  case_5 <- align_obj$cost_matrix[align_obj$T1_children[[i]][p], 'lambda'] + min(case_5)
  return(case_5)
} 

lambda_T22 <- function(align_obj, i, j) {
  subtree_cost <- c()
  T2_children_j <- align_obj$T2_children[[j]]
  for (r in T2_children_j) {
    subtree_cost <- c(subtree_cost, align_obj$T_matrix[i, r] - align_obj$T_matrix['lambda', r])
  }
  align_obj$T_matrix['lambda', j] + min(subtree_cost)
}

T1_lambda2 <- function(align_obj, i, j) {
  subtree_cost <- c()
  T1_children_i <- align_obj$T1_children[[i]]
  for (r in T1_children_i) {
    subtree_cost <- c(subtree_cost, align_obj$T_matrix[r, j] - align_obj$T_matrix[r, 'lambda'])
  }
  align_obj$T_matrix[i, 'lambda'] + min(subtree_cost)
}


T1_T22 <- function(align_obj, i, j) {
  n_i <- length(align_obj$T1_children[[i]])
  n_j <- length(align_obj$T2_children[[j]])
  F_i_j <- align_obj$F_map[[paste_collapse(1, i, 1, n_i, 2, j, 1, n_j)]] 
  F_i_j + align_obj$cost_matrix[i, j]
}




#' Traceback for constructing aligned tree
traceback2 <- function(align_obj) {
  align_obj$alignment <- c()
  align_obj <- recurse2(align_obj, align_obj$T1, align_obj$T2, left = T)
  return(align_obj)
}

recurse2 <- function(align_obj, x, y, left) {
  x_cond <- is.null(x)
  y_cond <- is.null(y)
  if (left & !(x_cond & y_cond)) align_obj$alignment <- c(align_obj$alignment, ')')
  
  if (x_cond & !y_cond) {
    align_obj$alignment <- postorder2(y, align_obj$alignment, left = F)
  } else if (y_cond & !x_cond) {
    align_obj$alignment <- postorder2(x, align_obj$alignment, left = T)
  } else if (!(x_cond | y_cond)) { 
    choix <- align_obj$T_choice[x$label, y$label]
    if (choix == 3) {
      align_obj$alignment <- c(align_obj$alignment,  paste_collapse(x$label, y$label))
      align_obj <- recurse(align_obj, x$left, y$left, left = T)
      align_obj <- recurse(align_obj, x$right, y$right, left = F)
    } else if (choix == 2) {
      align_obj$alignment <- c(align_obj$alignment, paste_collapse(x$label, 'lambda'))
      if (left) {
        align_obj <- recurse(align_obj, x$left, y, left = T)
        align_obj <- recurse(align_obj, x$right, NULL, left = F)
      } else {
        align_obj <- recurse(align_obj, x$left, NULL, left = T)
        align_obj <- recurse(align_obj, x$right, y, left = F)
      }
    } else if (choix == 1) {
      align_obj$alignment <- c(align_obj$alignment,  paste_collapse('lambda', y$label))
      if (left) {
        align_obj <- recurse(align_obj, x, y$left, left = T)
        align_obj <- recurse(align_obj, NULL, y$right, left = F)
      } else {
        align_obj <- recurse(align_obj, NULL, y$left, left = T)
        align_obj <- recurse(align_obj, x, y$right, left = F)
      }
    }
  }
  
  if (!left & !(x_cond & y_cond)) align_obj$alignment <- c(align_obj$alignment, '(')
  if (left & !(x_cond & y_cond)) align_obj$alignment <- c(align_obj$alignment, ',')
  return(align_obj)
}

serialize_helper2 <- function(x, vec, position) {
  n_children <- length(x$children)
  if (position == -1) vec <- c(vec, ')')
  vec <- c(vec, x$label)
  if (n_children == 0) return(vec)
  for (i in seq(n_children))  {
    position = ifelse(i == 1, -1, ifelse(i == n_children, 1, 0))
    vec <- serialize_helper2(x$children[[i]], vec, position = position)
    if (position < 1 & n_children > 1) vec <- c(vec, ',')
    if (position == 1) vec <- c(vec, '(')
  }
  return(vec)
}

postorder2 <- function(x, alignment, left) {
  if (is.null(x)) return(alignment)
  else {
    insert <- ifelse(left, paste_collapse(x$label, 'lambda'), paste_collapse('lambda', x$label))
    alignment <- c(alignment, insert)
    if (!is.null(x$left)) alignment <- c(alignment, ')')
    alignment <- postorder2(x$left, alignment, left)
    if (!is.null(x$right)) alignment <- c(alignment, ',')
    alignment <- postorder2(x$right, alignment, left)
    if (!is.null(x$left)) alignment <- c(alignment, '(')
  }
  return(alignment)
}



build_tree <- function(align_obj) {
  align_obj$alignment <- sapply(align_obj$alignment, function(x) gsub(' ', '_', x))
  text <- paste0(paste0(rev(align_obj$alignment[-c(1, length(align_obj$alignment))]), collapse = ''), ';')
  align_obj$tree <- as_binary_tree(ape::read.tree(text = text))
  return(align_obj)
}