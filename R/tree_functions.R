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

#' Convert a binary tree as phylo object
#' @param binary_tree A binary tree object
#' @return A phylo object
#' @export
as_phylo <- function(binary_tree) {
  text <- paste0(paste0(rev(serialize_helper(binary_tree, c(), T)[-1]), collapse = ''), ';')
  tree <- ape::read.tree(text = text)
  return(tree)
}


as_binary_tree <- function(x, ...) {
  UseMethod("as_binary_tree", x)
}

add_internal_node_label <- function(newick_string) {
  #print(newick_string)
  numbers <- na.omit(as.numeric(strsplit(newick_string, '[,();]+', fixed = F)[[1]]))
  i <- max(numbers) + 1
  vec <- strsplit(newick_string, '')[[1]]
  vec_ret <- c(vec[1])
  suppressWarnings(
    for (j in seq(2, length(vec))) {
      if (vec[j] %in% c(')', ';') & is.na(as.numeric(vec[j - 1]))) {
        vec_ret <- c(vec_ret, i)
        i <- i + 1
      }
      vec_ret <- c(vec_ret, vec[j])
    }
  )
  newick_string <- paste0(vec_ret, collapse = '')
  return(newick_string)
}

remove_branch_length <- function(newick_string) {
  gsub(':[0-9.]+', '', newick_string, fixed = F)
}

summarize_nodes <- function(tree, data, membership, f) {
  summary_stats_helper <- function(node, x) {
    
  }
  summary_stats <- c()
  summary_stats <- summary_stats_helper(tree, summary_stats)
  return(summary_stats)
}

#' Convert a hclust object as a binary tree object
#' @param hclust_obj A hclust object
#' @param data A data matrix
#' @param membership A membership matrix for leaves
#' @return A list containing the binary tree object and a vector of summary statistics on each node
#' @export
#' @import ape
as_binary_tree.hclust <- function(hclust_obj, data, membership, f) {
  tree <- ape::as.phylo(hclust_obj)
  newick_string <- add_internal_node_label(remove_branch_length(ape::write.tree(tree, digits = 0)))
  tree <- as_binary_tree.default(newick_string)
  summary_stats <- summarize_nodes(tree, data, membership, f)
  return(list(tree = tree, summary_stats = summary_stats))
} 

#' Convert a newick string as a binary tree object
#' @param newick_string A newick string
#' @return A binary tree object
#' @export
as_binary_tree.default <- function(newick_string) {
  vec <- convert_to_vec(newick_string)
  root <- binary_tree(label = vec[1])
  tree <- deserialize_helper(root, vec[-1])
  return(tree)
}

# convert phylo to binary tree object
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

# scan the vector to find matching right
break_vec <- function(vec) {
  # print(vec)
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

