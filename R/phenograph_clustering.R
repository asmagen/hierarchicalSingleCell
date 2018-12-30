#' Construct Hierarchy based on Cluster Mean
#'
#' Construct cluster hierarchy based on the mean of each cluster.
#'
#' @param data A data matrix to compute cluster mean..
#' @param membership A vector of cluster membership.
#' @return A hclust object.
#' @export
#' @author Meng Wang
construct_hierarchy <- function(data, membership, func) {
  cluster.mean = c()
  for (cluster in unique(membership)) {
    cluster.data = data[membership == cluster,]
    cluster.mean = rbind(cluster.mean, apply(cluster.data, 2, func))
  }
  rownames(cluster.mean) = as.character(unique(membership))
  dist = stats::dist(cluster.mean, method='euclidean')
  hclust = stats::hclust(dist, method='complete')
  dend = as.dendrogram(hclust)
  dend_k = dendextend::find_k(dend)
  dend = dendextend::color_branches(dend, dend_k$k)
  dend2 = dendextend::seriate_dendrogram(dend, dist)
  ret = list(data=cluster.mean, dend=dend2)
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
