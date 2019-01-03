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

#' Generic function for plotting binary tree
plot <- function(x, ...) {
  UseMethod("plot", x)
}

#' Plot binary tree
#' @param tree A binary tree object
#' @export
plot.binary_tree <- function(tree) {
  plot(as_phylo(tree), show.node.label = T, no.margin = F, 
       edge.width = 2, edge.color = 'gray', cex = 1.4, font = 1, 
       label.offset = 0.1)
}

#' Align two binary ordered trees
#' @param T1 A binary tree
#' @param T2 A binary tree
#' @param cost_matrix Cost matrix between each pair of nodes in T1 and T2
#' @return An alignment object
#' @export
align <- function(T1, T2, cost_matrix) {
  align_obj <- create_align_object(T1, T2, cost_matrix)
  align_obj <- initialize(align_obj)
  align_obj <- fill_matrix(align_obj)
  align_obj <- traceback(align_obj)
  align_obj <- build_tree(align_obj)
  return(align_obj)
}

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
  T1_node <- rownames(align_object$T_matrix)[-1]
  for (node in T1_node) {
    children <- align_object$T1_children[node, ]
    F_node_lambda <- align_object$T_matrix[children, "lambda"]
    align_object$F_map[[paste_collapse(1, node, 1, 1, 'lambda')]] <- 2 * F_node_lambda[1]
    align_object$F_map[[paste_collapse(1, node, 2, 2, 'lambda')]] <- 2 * F_node_lambda[2]
    align_object$F_map[[paste_collapse(1, node, 1, 2, 'lambda')]] <- sum(F_node_lambda)
    align_object$T_matrix[node, "lambda"] <- sum(F_node_lambda) + align_object$cost_matrix[node, "lambda"]
  }
  return(align_object)
}

initialize_T2 <- function(align_object) {
  T2_node <- colnames(align_object$T_matrix)[-1]
  for (node in T2_node) {
    children <- align_object$T2_children[node, ]
    lambda_F_node <- align_object$T_matrix["lambda", children]
    align_object$F_map[[paste_collapse('lambda', 2, node, 1, 1)]] <- 2 * lambda_F_node[1]
    align_object$F_map[[paste_collapse('lambda', 2, node, 2, 2)]] <- 2 * lambda_F_node[2]
    align_object$F_map[[paste_collapse('lambda', 2, node, 1, 2)]] <- sum(lambda_F_node)
    align_object$T_matrix["lambda", node] <- sum(lambda_F_node) + align_object$cost_matrix["lambda", node]
  }
  return(align_object)
}

#' Intialize the T matrix and F map
#' 
#' @param align_object An alignment object
#' @return An alignment object with T matrix and F map initialized
#' @export
initialize <- function(align_object) {
  align_object$T1_children <- create_children_map(align_object$T1)
  align_object$T2_children <- create_children_map(align_object$T2)
  align_object <- initialize_T1(align_object)
  align_object <- initialize_T2(align_object)
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

paste_collapse <- function(...) paste0(c(...), collapse = '_')


T1_T2 <- function(align_obj, i, j) {
  F_i_j <- align_obj$F_map[[paste_collapse(1, i, 1, 2, 2, j, 1, 2)]] 
  F_i_j + align_obj$cost_matrix[i, j]
}

#' Fill T matrix and F map
#' 
#' @param align_obj An alignment object
#' @return An alignment object with T matrix and F map filled
#' @export
fill_matrix <- function(align_obj) {
  align_obj$T_choice <- align_obj$T_matrix
  align_obj$T_choice[,] <- NA
  
  for (i in rownames(align_obj$T_matrix)[-1]) {
    for (j in colnames(align_obj$T_matrix)[-1]) {
      cost <- c()
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
      # print(cost)
      align_obj$T_choice[i, j] <- which.min(cost)
      align_obj$T_matrix[i, j] <- min(cost, na.rm = T)
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

fill_case_4 <- function(align_obj, i, s, p, j, t, q) {
  case_4 <- c()
  for (k in seq(s, p)) {
    idx1 <- paste_collapse(1, i, k, p)
    T2_j_q <- align_obj$T2_children[j, q]
    if (T2_j_q == 'lambda') idx2 <- 'lambda'
    else idx2 <- paste_collapse(2, T2_j_q, 1, 2)
    case_4 <- c(case_4, align_obj$F_map[[paste_collapse(1, i, s, k - 1, 2, j, t, q - 1)]] + align_obj$F_map[[paste_collapse(idx1, idx2)]])
  }
  case_4 <- align_obj$cost_matrix['lambda', align_obj$T2_children[j, q]] + min(case_4)
  return(case_4)
}

fill_case_5 <- function(align_obj, i, s, p, j, t, q) {
  case_5 <- c()
  for (k in seq(t, q)) {
    T1_i_p <- align_obj$T1_children[i, p]
    if (T1_i_p == 'lambda') idx1 <- 'lambda'
    else idx1 <- paste_collapse(1, T1_i_p, 1, 2)
    idx2 <- paste_collapse(2, j, k, q)
    case_5 <- c(case_5, align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, k - 1)]] + align_obj$F_map[[paste_collapse(idx1, idx2)]])
  }
  case_5 <- align_obj$cost_matrix[align_obj$T1_children[i, p], 'lambda'] + min(case_5)
  return(case_5)
} 

fill_F_helper <- function(align_obj, i, s, p, j, t, q) {
  case_1 <- align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, q)]] +
    align_obj$T_matrix[align_obj$T1_children[i, p], 'lambda']
  case_2 <- align_obj$F_map[[paste_collapse(1, i, s, p, 2, j, t, q - 1)]] +
    align_obj$T_matrix['lambda', align_obj$T2_children[j, q]]
  case_3 <- align_obj$F_map[[paste_collapse(1, i, s, p - 1, 2, j, t, q - 1)]] +
    align_obj$T_matrix[align_obj$T1_children[i, p], align_obj$T2_children[j, q]]
  case_4 <- fill_case_4(align_obj, i, s, p, j, t, q)
  case_5 <- fill_case_5(align_obj, i, s, p, j, t, q)

  all <- c(case_1, case_2, case_3, case_4, case_5)
  loc <- paste_collapse(1, i, s, p, 2, j, t, q)
  align_obj$F_map[[loc]] <- min(all)
  return(align_obj)
}


#' Traceback for constructing aligned tree
traceback <- function(align_obj) {
  align_obj$alignment <- c()
  align_obj <- recurse(align_obj, align_obj$T1, align_obj$T2, left = T)
  return(align_obj)
}

recurse <- function(align_obj, x, y, left) {
  x_cond <- is.null(x)
  y_cond <- is.null(y)
  if (left & !(x_cond & y_cond)) align_obj$alignment <- c(align_obj$alignment, ')')
  if (x_cond & !y_cond) {
    postorder(y, function(t) {
      align_obj$alignment <<- c(align_obj$alignment, paste_collapse('lambda', t$label))
    })
  } else if (y_cond & !x_cond) {
    postorder(x, function(t) {
      align_obj$alignment <<- c(align_obj$alignment, paste_collapse(t$label, 'lambda'))
    })
  } 
  else if (!(x_cond | y_cond)) {
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

build_tree <- function(align_obj) {
  text <- paste0(paste0(rev(align_obj$alignment[-c(1, length(align_obj$alignment))]), collapse = ''), ';')
  align_obj$tree <- ape::read.tree(text = text)
  return(align_obj)
}

#' convert phylo to binary tree object
serialize_helper <- function(x, vec, left) {
  if (is.null(x)) return(vec)
  else {
    if (left) vec <- c(vec, ')')
    vec <- c(vec, x$label)
    vec <- serialize_helper(x$left, vec, left = T)
    if (!is.null(x$left)) vec <- c(vec, ',')
    vec <- serialize_helper(x$right, vec, left = F)
    if (!is.null(x$right)) vec <- c(vec, '(')
  }
  return(vec)
}

as_phylo <- function(binary_tree) {
  text <- paste0(paste0(rev(serialize_helper(binary_tree, c(), T)[-1]), collapse = ''), ';')
  tree <- ape::read.tree(text = text)
  return(tree)
}

#' scan the vector to find matching right
break_vec <- function(vec) {
  print(vec)
  found <- 0
  if (vec[1] == ',' && length(vec) == 3) return(list(NULL, vec[2]))
  for (i in seq_along(vec)) {
    if (vec[i] == ',' && found == 0) {
      left <- vec[1:(i-1)]
      right <- vec[(i+1):(length(vec) - 1)]
      return(list(left, right))
    } else {
      if (vec[i] == ')') found = found + 1
      if (vec[i] == '(') found = found - 1
    }
  }
}

deserialize_helper <- function(x, vec) {
  if (length(vec) == 0) return(x)
  else if (vec[1] == ')') {
    broken_vec <- break_vec(vec[-seq(2)])
    x$left <- deserialize_helper(binary_tree(label = vec[2]), broken_vec[[1]])
    x$right <- deserialize_helper(binary_tree(label = broken_vec[[2]][1]), broken_vec[[2]][-1])
  }
  return(x)
}

as_binary_tree <- function(newick_string) {
  vec <- convert_to_vec(newick_string)
  root <- binary_tree(label = vec[1])
  tree <- deserialize_helper(root, vec[-1])
  return(tree)
}

convert_to_vec <- function(newick_string) {
  vec <- c()
  all_char <- strsplit(newick_string, split = '')[[1]]
  temp <- c()
  for (i in seq_along(all_char)) {
    if (all_char[i] %in% c(',', '(', ')', ';')) {
      vec <- c(vec, paste0(temp, collapse = ''), all_char[i])
      temp <- c()
    }
    else temp <- c(temp, all_char[i])
  }
  vec <- rev(vec)[-c(1, length(vec))]
  vec <- vec[vec != '']
  return(vec)
}
