#' Construct Hierarchy based on Cluster Mean
#'
#' Construct cluster hierarchy based on the mean of each cluster.
#'
#' @param data A data matrix to compute cluster mean..
#' @param membership A vector of cluster membership.
#' @return A hclust object.
#' @export
#' @import purrr
construct_hierarchy <- function(data, membership, dist_method = 'euclidean', method = 'mean', 
                                hclust_method = 'complete', foldchange = F) {
  func <- switch(method,
                 mean   = function(v) log(base::mean(exp(v) - 1) + 1),
                 median = function(v) log(stats::median(exp(v) - 1) + 1))
  
  cluster_summary <- c()
  if (foldchange) {
    clusters <- unique(membership)
    for (i in clusters) {
      cluster_summary <- rbind(cluster_summary, as.numeric(map_dbl(data[membership == i,], func) / map_dbl(data[membership != i,], func)))
    }
  } else {
    cluster_split <- split(data.frame(data), membership)
    clusters <- names(cluster_split)
    cluster_split_summary <- map(cluster_split, ~ map_dbl(.x, func))
    walk(cluster_split_summary, function(x) cluster_summary <<- rbind(cluster_summary, x))
  }
  rownames(cluster_summary) <- clusters
  
  dist_func <- switch(dist_method,
                      euclidean = function(x) stats::dist(x, method='euclidean'),
                      spearman  = function(x) as.dist(1 - cor(t(x), method = 'spearman')),
                      pearson   = function(x) as.dist(1 - cor(t(x), method = 'pearson')))
  
  dist <- dist_func(cluster_summary)
  hclust <- stats::hclust(dist, method = hclust_method)
  dend <- as.dendrogram(hclust)
  dend_k <- dendextend::find_k(dend)
  dend <- dendextend::color_branches(dend, dend_k$k)
  dend2 <- dendextend::seriate_dendrogram(dend, dist)
  ret <- list(dend=dend2, hclust = hclust, data = data, membership = membership)
  class(ret) <- "hierarchy"
  return(ret)
}

#' Plot Hierarchy Constructed
#'
#' @param x A hierarchy object.
#' @return A heatmap.2 object
#' @export
plot.hierarchy = function(x, ha_row=NULL, show_row_names=T, show_col_names=F) {
  d.marker = stats::dist(t(x$data), method='euclidean')
  hc.marker = stats::hclust(d.marker, method='ward.D2')
  dend = stats::as.dendrogram(hc.marker)
  dend_k = dendextend::find_k(dend)
  dend = dendextend::color_branches(dend, dend_k$k)
  dend2 = dendextend::seriate_dendrogram(dend, d.marker)
  if (!is.null(ha_row)) {
    ComplexHeatmap::Heatmap(x$data, cluster_columns=dend2, cluster_rows=x$dend,
                            col=viridis::viridis(200), show_row_names=show_row_names,
                            show_column_names = show_col_names) + ha_row
  } else {
    ComplexHeatmap::Heatmap(x$data, cluster_columns=dend2, cluster_rows=x$dend,
                            col=viridis::viridis(200), show_row_names = show_row_names,
                            show_column_names = show_col_names)
  }
  
}

#' Construct clustering tree from proportion graph
#' @param clusters.aggregate A data matrix with cluster membership at different levels (column names should be K1, K2, ...)
#' @return A list containing a proportion graph and a trimmed version keeping only the incident edge with highest weight for each node
#' @export
#' @import tidygraph
clustertree_construction <- function(clusters.aggregate) {
  prefix = "K"
  node_aes_list <- list(colour = list(value = prefix, aggr = NULL), size = list(value = "size", aggr = NULL), 
                        alpha = list(value = 1, aggr = NULL))
  graph <- clustree:::build_tree_graph(clusters.aggregate, "K", metadata = NULL, count_filter = 0, prop_filter = 0.1, node_aes_list = node_aes_list)
  graph_trimmed <- graph %>%
    tidygraph::activate(edges) %>% 
    tidygraph::group_by(to_K, to_clust) %>%
    tidygraph::top_n(1, in_prop)
  return(list(graph = graph, trimmed = graph_trimmed))
}

#' Plot a igraph object
#' @param graph An igraph object
#' @export
#' @import ggraph
plot_graph <- function(graph) {
  ggraph::ggraph(graph, 'igraph', algorithm = 'tree') + 
    geom_edge_link() +
    ggforce::theme_no_axes()
}