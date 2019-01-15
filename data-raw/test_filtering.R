library(ggplot2,quietly=T)
data('T3K')
data('T4K')

# examine the proportion of ribosomal proteins in the dataset
# defined as starting with RP followed by numbers
# 10X : isolated Pan T cells 40 - 45% of reads map to RP (https://kb.10xgenomics.com/hc/en-us/articles/218169723-What-fraction-of-reads-map-to-ribosomal-proteins-)
# startegy 1: only include variable genes 
# check the amount of ribsosomal proteins in the dataset
dir_path <- '~/Desktop'
HVG_names <- ReSET::get_HVG(T3K$normalized_data, dir_path, min.mean = 0.3, min.frac.cells = 5, min.dispersion.scaled = 1,
                            remove_RP = T)
filtered_T3K <- list(normalized_data = T3K$normalized_data[HVG_names,], membership = T3K$membership)

# repeat the same procedure for T4K dataset
HVG_names <- ReSET::get_HVG(T4K$normalized_data, dir_path, min.mean = 0.3, min.frac.cells = 5, min.dispersion.scaled = 1,
                            remove_RP = T)
filtered_T4K <- list(normalized_data = T4K$normalized_data[HVG_names,], membership = T4K$membership)


configurations <- list(c('euclidean', 'complete'), c('spearman', 'complete'), c('pearson', 'complete'),
                       c('euclidean', 'average'), c('spearman', 'average'), c('pearson', 'average'))
lapply(configurations, function(configuration) {
  png(file = paste0('~/Desktop/', paste0(configuration, collapse = '_'), '.png'))
  par(mfrow= c(1, 2)); ReSET::test_hierarchy_construction(filtered_T3K, filtered_T4K, configuration, foldchange = T)
  dev.off()
})
lapply(configurations, function(configuration) {
  png(file = paste0('~/Desktop/', paste0(configuration, collapse = '_'), '.png'))
  par(mfrow= c(1, 2)); ReSET::test_hierarchy_construction(filtered_T3K, filtered_T4K, configuration, foldchange = F)
  dev.off()
})

T3K_tree <- ReSET::construct_hierarchy(t(filtered_T3K$normalized_data), filtered_T3K$membership, dist_method = 'pearson', foldchange = F)
T4K_tree <- ReSET::construct_hierarchy(t(filtered_T4K$normalized_data), filtered_T4K$membership, dist_method = 'pearson', foldchange = F)

# The conclusion is that using highly variable genes help with hierarchy construction
# The most suitable hierarchy is dependent on the metrics and data at hand 
# for now the best model for T3K is (mean, pearson/pearman) which exactly where clusters 
# for the same cell type are closest to each other
# for T4K is harder to interpret

# Let use the two constructed hierarchy to perform the alignment
summary_method <-  'mean'
T3K_binary_tree <- with(T3K_tree, ReSET::as_binary_tree(hclust, data, membership, summary_method))
T4K_binary_tree <- with(T4K_tree, ReSET::as_binary_tree(hclust, data, membership, summary_method))
cost_matrix1 <- ReSET::cor_cost_matrix(T3K_binary_tree$summary_stats, T4K_binary_tree$summary_stats, method = 'spearman')
cost_matrix2 <- ReSET::cor_cost_matrix(T3K_binary_tree$summary_stats, T4K_binary_tree$summary_stats, method = 'pearson')
cost_matrix3 <- ReSET::euc_cost_matrix(T3K_binary_tree$summary_stats, T4K_binary_tree$summary_stats)

library(iheatmapr)
ReSET:::plot_heatmap(1 - cost_matrix1, name = 'spearman')
ReSET:::plot_heatmap(1 - cost_matrix2, name = 'pearson')
ReSET:::plot_heatmap(-cost_matrix3, name = '-euclidean')

aligned_tree1 <- ReSET::align(T3K_binary_tree$tree, T4K_binary_tree$tree, cost_matrix1)
aligned_tree2 <- ReSET::align(T3K_binary_tree$tree, T4K_binary_tree$tree, cost_matrix2)
aligned_tree3 <- ReSET::align(T3K_binary_tree$tree, T4K_binary_tree$tree, cost_matrix3)

options(mfrow = c(1, 3))
plot(aligned_tree1$tree)
plot(aligned_tree2$tree)
plot(aligned_tree3$tree)

# conclusion: if the hierarchy constructed is not matching, then the alignment will be hard
# 1. tune the parameters so that the trees resemble each other  (find a better method 
# to construct hierarchy of trees)
# 2. what are some hierarchy construction method (tuning on the distance matrix)
# distance in the 
# 3. fold change doesn't seem to change the resolution 