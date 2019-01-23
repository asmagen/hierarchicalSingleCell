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
#' @export
plot <- function(x, ...) {
  UseMethod("plot", x)
}

#' Plot binary tree
#' @param tree A binary tree object
#' @export
plot.binary_tree <- function(tree, ...) {
  ape::plot.phylo(as_phylo(tree), show.node.label = T, no.margin = F, 
       edge.width = 2.3, edge.color = 'gray', cex = 0.9, font = 1, 
       label.offset = 0.1, underscore = T, ...)
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


#' Generic method for convert trees to binary tree
#' 
#' `as_binary_tree` is a generic function for converting various objects into a 
#' binary_tree object
#' @param object An object to be converted to a binary_tree
#' @param ... Additional arguments for conversion
#' @export
as_binary_tree <- function(x, ...) {
  UseMethod("as_binary_tree", x)
}

add_internal_node_label <- function(newick_string) {
  labels <- strsplit(newick_string, '[,();]+', fixed = F)[[1]]
  labels <- labels[labels != '']
  suppressWarnings(numbers <- na.omit(as.numeric(labels)))
  i <- ifelse(length(numbers)!=0, max(numbers) + 1, 1)
  vec <- strsplit(newick_string, '')[[1]]
  vec_ret <- c(vec[1])
  suppressWarnings(
    for (j in seq(2, length(vec))) {
      if (vec[j] %in% c(')', ';', ',') & (vec[j - 1] == ')')) {
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
  gsub(':[0-9.e+]+', '', newick_string, fixed = F)
}

get_all_children <- function(x) {
  if (is.null(x$left) & is.null(x$right)) {
    return(x$label)
  } else {
    return(c(get_all_children(x$left), get_all_children(x$right)))
  }
}

get_all_children_leaves <- function(tree, children_map) {
  children_map[[tree$label]] <- get_all_children(tree)
  if (!is.null(tree$left)) children_map <- get_all_children_leaves(tree$left, children_map)
  if (!is.null(tree$right)) children_map <- get_all_children_leaves(tree$right, children_map)
  return(children_map)
}

summarize_nodes <- function(tree, data, membership, f) {
  children_map <- get_all_children_leaves(tree, list())
  
  summary_stats_helper <- function(node, x) {
    if(is.null(node)) return(x)
    summarized <- apply(data[membership %in% children_map[[node$label]], ], 2, f)
    x <- rbind(x, c(node$label, summarized))
    x <- summary_stats_helper(node$left, x)
    x <- summary_stats_helper(node$right, x)
    return(x)
  }
  summary_stats_helper(tree, c())
}

#' Convert a hclust object as a binary tree object
#' @param hclust_obj A hclust object
#' @param data A data matrix
#' @param membership A membership matrix for leaves
#' @param method A string indicating method for summarizing the nodes
#' @return A list containing the binary tree object and a vector of summary statistics on each node
#' @export
#' @import ape
as_binary_tree.hclust <- function(hclust_obj, data, membership, method = 'median') {
  f <- switch(method,
              median = function(v) log(stats::median(exp(v) - 1) + 1),
              mean = function(v) log(base::mean(exp(v) - 1) + 1))
  tree <- ape::as.phylo(hclust_obj)
  newick_string <- add_internal_node_label(remove_branch_length(ape::write.tree(tree, digits = 0)))
  tree <- as_binary_tree.default(newick_string)
  summary_stats <- summarize_nodes(tree, data, membership, f)
  cluster <- unname(summary_stats[,1])
  summary_stats <- summary_stats[,-1]
  rownames(summary_stats) <- cluster
  class(summary_stats) <- 'numeric'
  return(list(tree = tree, summary_stats = summary_stats))
}

#' Convert a phylo object as a binary tree object
#' @param phylo_obj A phylo object
#' @return A binary tree object 
#' @export
as_binary_tree.phylo <- function(phylo_obj) {
  newick_string <- ape::write.tree(phylo_obj, digits = 0)
  tree <- as_binary_tree.default(newick_string)
  return(tree)
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

