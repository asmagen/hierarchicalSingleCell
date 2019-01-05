library(dplyr)

# preprocess the dataset
T3K_path <- "inst/extdata/3k_pan_T_filtered_gene_bc_matrices/GRCh38"
T4K_path <- "inst/extdata/4k_pan_T_filtered_gene_bc_matrices/GRCh38"

T3K_data <- ReSET::preprocess_data(T3K_path)
T4K_data <- ReSET::preprocess_data(T4K_path)

# construct hierarchy
T3K_data <- Seurat::FindClusters(T3K_data, reduction.type = "pca", dims.use = 1:12, resolution = 2, print.output = 0) # 12 clusters
T3K_hierarchy <- ReSET::construct_hierarchy(t(T3K_data@scale.data[T3K_data@var.genes, ]), membership = T3K_data@ident, mean)
stats:::plot.dendrogram(T3K_hierarchy$dend)

T4K_data <- Seurat::FindClusters(T4K_data, reduction.type = "pca", dims.use = 1:12, resolution = 2, print.output = 0) # 13 clusters
T4K_hierarchy <- ReSET::construct_hierarchy(t(T4K_data@scale.data[T4K_data@var.genes, ]), membership = T4K_data@ident, mean)
stats:::plot.dendrogram(T4K_hierarchy$dend)

# Let's try it on real dataset
T3K_tree <- with(T3K_hierarchy, ReSET::as_binary_tree(hclust, data, membership, mean))
T4K_tree <- with(T4K_hierarchy, ReSET::as_binary_tree(hclust, data, membership, mean))

cost_matrix <- matrix(NA, nrow = nrow(T3K_tree$summary_stats) + 1, ncol = nrow(T4K_tree$summary_stats) + 1)
rownames(cost_matrix) <- c('lambda', rownames(T3K_tree$summary_stats))
colnames(cost_matrix) <- c('lambda', rownames(T4K_tree$summary_stats))
intersect_cols <- intersect(colnames(T3K_tree$summary_stats), colnames(T4K_tree$summary_stats))
euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

for (i in rownames(cost_matrix)) {
  for (j in colnames(cost_matrix)) {
    if (i == 'lambda' & j == 'lambda') cost_matrix[i, j] <- 0
    else if (i == 'lambda' | j == 'lambda') cost_matrix[i, j] <- 2
    else cost_matrix[i, j] <- euc_dist(T3K_tree$summary_stats[i,intersect_cols], T4K_tree$summary_stats[j,intersect_cols])
  }
}



#aligned <- ReSET::align(T3K_tree$tree, T4K_tree$tree, cost_matrix)