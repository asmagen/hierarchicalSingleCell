## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
T3K_path <- system.file("extdata", '3k_pan_T_filtered_gene_bc_matrices/GRCh38', package = 'ReSET', mustWork = T)
T4K_path <- system.file("extdata", "4k_pan_T_filtered_gene_bc_matrices/GRCh38", package = 'ReSET', mustWork = T)

T3K_data <- ReSET::preprocess_data(T3K_path)
T4K_data <- ReSET::preprocess_data(T4K_path)

## ----message=FALSE, warning=FALSE, fig.width = 6, fig.height = 3---------
T3K_data <- Seurat::FindClusters(T3K_data, reduction.type = "pca", dims.use = 1:12, resolution = 2, print.output = 0) # 12 clusters
T3K_hierarchy <- ReSET::construct_hierarchy(T3K_data@scale.data, membership = T3K_data@ident, mean)
plot(T3K_hierarchy$dend)

## ----message=FALSE, warning=FALSE, fig.width = 6, fig.height = 3---------
T4K_data <- Seurat::FindClusters(T4K_data, reduction.type = "pca", dims.use = 1:12, resolution = 2, print.output = 0) # 13 clusters
T4K_hierarchy <- ReSET::construct_hierarchy(T4K_data@scale.data, membership = T4K_data@ident, mean)
plot(T4K_hierarchy$dend)

