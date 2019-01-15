#' Test hierarchy construction
#' @param T3K a hierarchy object (for example T3K)
#' @param T4K A hierarchy object (for example T4K)
#' @param configuration A list of configuration for hierarchy construction
#' @export
test_hierarchy_construction <- function(T3K, T4K, configuration, foldchange = T) {
  #print(configuration)
  T3K_tree <- ReSET::construct_hierarchy(t(T3K$normalized_data), T3K$membership, dist_method = configuration[1],
                                         hclust_method = configuration[2], foldchange = foldchange)
  T4K_tree <- ReSET::construct_hierarchy(t(T4K$normalized_data), T4K$membership, dist_method = configuration[1],
                                         hclust_method = configuration[2], foldchange = foldchange)
  stats:::plot.dendrogram(T3K_tree$dend, main = paste0(configuration, collapse = ', '))
  stats:::plot.dendrogram(T4K_tree$dend, main = paste0(configuration, collapse = ', '))
  return(list(T3K_tree, T4K_tree))
}

#' Test alignment method
#' @param T3K a hierarchy object (for example T3K)
#' @param T4K A hierarchy object (for example T4K)
#' @param configuration A list of configuration for hierarchy construction
#' @export
test_alignment_configuration <- function(T3K_tree, T4K_tree, configuration) {
  T3K_binary_tree <- with(T3K_tree, ReSET::as_binary_tree(hclust, data, membership, configuration[1]))
  T4K_binary_tree <- with(T4K_tree, ReSET::as_binary_tree(hclust, data, membership, configuration[1]))
  cost_matrix <- switch(configuration[2],
                        euclidean = ReSET::euc_cost_matrix(T3K_binary_tree$summary_stats, T4K_binary_tree$summary_stats),
                        pearson = ReSET::cor_cost_matrix(T3K_binary_tree$summary_stats, T4K_binary_tree$summary_stats, method = 'pearson'),
                        spearman = ReSET::cor_cost_matrix(T3K_binary_tree$summary_stats, T4K_binary_tree$summary_stats, method = 'spearman'))
  align_obj <- ReSET::align(T3K_binary_tree$tree, T4K_binary_tree$tree, cost_matrix)
  plot(align_obj$tree, main = paste0(configuration, collapse = ', '))
}

plot_heatmap <- function(dist_mat, name) {
  main_heatmap(dist_mat[-1,-1], name = name) %>%
    add_row_labels() %>%
    add_col_labels() %>%
    add_row_clustering() %>%
    add_col_clustering()
}