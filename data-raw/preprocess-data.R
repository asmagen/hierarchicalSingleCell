library(dplyr)

# preprocess the dataset
T3K_path <- "inst/extdata/3k_pan_T_filtered_gene_bc_matrices/GRCh38"
T4K_path <- "inst/extdata/4k_pan_T_filtered_gene_bc_matrices/GRCh38"

T3K_data <- ReSET::preprocess_data(T3K_path)
T4K_data <- ReSET::preprocess_data(T4K_path)

# construct hierarchy
T3K_data <- Seurat::FindClusters(T3K_data, reduction.type = "pca", dims.use = 1:12, resolution = 1, print.output = 0) 
T3K_hierarchy <- ReSET::construct_hierarchy(t(T3K_data@scale.data[T3K_data@var.genes, ]), membership = T3K_data@ident, mean)
stats:::plot.dendrogram(T3K_hierarchy$dend)

T4K_data <- Seurat::FindClusters(T4K_data, reduction.type = "pca", dims.use = 1:12, resolution = 1, print.output = 0) # 13 clusters
T4K_hierarchy <- ReSET::construct_hierarchy(t(T4K_data@scale.data[T4K_data@var.genes, ]), membership = T4K_data@ident, mean)
stats:::plot.dendrogram(T4K_hierarchy$dend)

# Let's try it on real dataset
T3K_tree <- with(T3K_hierarchy, ReSET::as_binary_tree(hclust, data, membership, mean))
T4K_tree <- with(T4K_hierarchy, ReSET::as_binary_tree(hclust, data, membership, mean))

par(mfrow = c(1, 2))
plot(T3K_tree$tree)
plot(T4K_tree$tree)

cost_matrix <- matrix(NA, nrow = nrow(T3K_tree$summary_stats) + 1, ncol = nrow(T4K_tree$summary_stats) + 1)
rownames(cost_matrix) <- c('lambda', rownames(T3K_tree$summary_stats))
colnames(cost_matrix) <- c('lambda', rownames(T4K_tree$summary_stats))
intersect_cols <- intersect(colnames(T3K_tree$summary_stats), colnames(T4K_tree$summary_stats))
euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

for (i in rownames(cost_matrix)) {
  for (j in colnames(cost_matrix)) {
    if (i == 'lambda' & j == 'lambda') cost_matrix[i, j] <- 0
    else if (i == 'lambda' | j == 'lambda') cost_matrix[i, j] <- 1
    else cost_matrix[i, j] <- 1 - cor(T3K_tree$summary_stats[i,intersect_cols], T4K_tree$summary_stats[j,intersect_cols])
  }
}

# align_obj <- create_align_object(T3K_tree$tree, T4K_tree$tree, cost_matrix) # align with itself
# align_obj <- initialize(align_obj)
# align_obj <- fill_matrix(align_obj)
# align_obj <- traceback(align_obj)
# align_obj <- build_tree(align_obj)
align_obj <- ReSET::align(T3K_tree$tree, T4K_tree$tree, cost_matrix)
plot(as_binary_tree(align_obj$tree))


# Let's try the annotated dataset from Assaf
load('inst/extdata/T3k.RData')

T3K_tree2 <- ReSET::construct_hierarchy(t(normalized), tree[,9], median)
labels <- unique(tree[,c(8,9)])
T3K_map <- hashmap::hashmap(labels[,1], labels[,2])

load('inst/extdata/T4k.RData')
T4K_tree2 <- ReSET::construct_hierarchy(t(normalized), tree[,11], median)
labels <- unique(tree[,c(10,11)])
T4K_map <- hashmap::hashmap(labels[,1], labels[,2])

stats:::plot.dendrogram(T3K_tree2$dend)
stats:::plot.dendrogram(T4K_tree2$dend)

T3K_binary_tree <- with(T3K_tree2, ReSET::as_binary_tree(hclust, data, membership, median))
T4K_binary_tree <- with(T4K_tree2, ReSET::as_binary_tree(hclust, data, membership, median))
par(mfrow = c(1, 2))
plot(T3K_binary_tree$tree)
plot(T4K_binary_tree$tree)

cost_matrix <- ReSET::cor_cost_matrix(T3K_binary_tree$summary_stats, T4K_binary_tree$summary_stats)
cost_matrix[1:5, 1:5]
align_obj <- ReSET::align(T3K_binary_tree$tree, T4K_binary_tree$tree, cost_matrix)
plot(align_obj$tree)

library(iheatmapr)
main_heatmap(1 - cost_matrix[-1,-1], name = 'correlation') %>%
  add_row_labels() %>%
  add_col_labels() %>%
  add_row_clustering() %>%
  add_col_clustering()
