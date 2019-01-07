#' Construct Hierarchy based on Cluster Mean
#'
#' Construct cluster hierarchy based on the mean of each cluster.
#'
#' @param data A data matrix to compute cluster mean..
#' @param membership A vector of cluster membership.
#' @return A hclust object.
#' @export
construct_hierarchy <- function(data, membership, method = 'median', dist_method = 'euclidean',
                                hclust_method = 'complete') {
  if (method == 'mean') func <- function(v) log(base::mean(exp(v) - 1) + 1)
  else if (method == 'median') func <- function(v) log(stats::median(exp(v) - 1) + 1)
  cluster_split <- lapply(split(data.frame(data), membership), function(x) apply(x, 2, func))
  clusters <- names(cluster_split)
  cluster_summary <- c()
  dump <- lapply(cluster_split, function(x) cluster_summary <<- rbind(cluster_summary, x))
  rownames(cluster_summary) <- clusters
  if (dist_method == 'euclidean') {
    dist <- stats::dist(cluster_summary, method='euclidean')
  } else if (dist_method == 'spearman') {
    dist <- as.dist(1 - cor(t(cluster_summary), method = 'spearman'))
  } else if (dist_method == 'pearson') {
    dist <- as.dist(1 - cor(t(cluster_summary), method = 'pearson'))
  }
  hclust = stats::hclust(dist, method = hclust_method)
  dend = as.dendrogram(hclust)
  dend_k = dendextend::find_k(dend)
  dend = dendextend::color_branches(dend, dend_k$k)
  dend2 = dendextend::seriate_dendrogram(dend, dist)
  ret = list(dend=dend2, hclust = hclust, data = data, membership = membership)
  class(ret) = "hierarchy"
  return(ret)
}

#' Plot Hierarchy Constructed
#'
#' @param x A hierarchy object.
#' @return A heatmap.2 object
#' @export
#' @author Meng Wang
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
