## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
T3K_path <- system.file("extdata", '3k_pan_T_filtered_gene_bc_matrices/GRCh38', package = 'ReSET', mustWork = T)
T4K_path <- system.file("extdata", "4k_pan_T_filtered_gene_bc_matrices/GRCh38", package = 'ReSET', mustWork = T)

T3K_data <- ReSET::preprocess_data(T3K_path)
T4K_data <- ReSET::preprocess_data(T4K_path)

## ----message=FALSE, warning=FALSE, fig.width = 6, fig.height = 3---------
T3K_data <- Seurat::FindClusters(T3K_data, reduction.type = "pca", dims.use = 1:12, resolution = 2, print.output = 0) # 12 clusters
T3K_hierarchy <- ReSET::construct_hierarchy(t(T3K_data@scale.data[T3K_data@var.genes, ]), membership = T3K_data@ident, mean)
stats:::plot.dendrogram(T3K_hierarchy$dend)

## ----message=FALSE, warning=FALSE, fig.width = 6, fig.height = 3---------
T4K_data <- Seurat::FindClusters(T4K_data, reduction.type = "pca", dims.use = 1:12, resolution = 2, print.output = 0) # 13 clusters
T4K_hierarchy <- ReSET::construct_hierarchy(t(T4K_data@scale.data[T4K_data@var.genes, ]), membership = T4K_data@ident, mean)
stats:::plot.dendrogram(T4K_hierarchy$dend)

## ----message=FALSE, warning=FALSE----------------------------------------
T3K_tree <- with(T3K_hierarchy, ReSET::as_binary_tree(hclust, data, membership, mean))
T4K_tree <- with(T4K_hierarchy, ReSET::as_binary_tree(hclust, data, membership, mean))

## ----message=FALSE, warning=FALSE, fig.width = 8, fig.height = 4---------
par(mfrow = c(1, 2))
plot(T3K_tree$tree)
plot(T4K_tree$tree)

## ----message=FALSE, warning=FALSE, fig.width = 8, fig.height = 4---------
cost_matrix <- ReSET::cor_cost_matrix(T3K_tree$summary_stats, T4K_tree$summary_stats)
cost_matrix[1:5, 1:5]

## ----message=FALSE, warning=FALSE, fig.width = 8, fig.height = 6---------
align_obj <- ReSET::align(T3K_tree$tree, T4K_tree$tree, cost_matrix)
plot(align_obj$tree)

