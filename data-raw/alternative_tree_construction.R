# use clustree to construct the hierarchy
library(clustree)
library(Rphenograph)
library(ggraph)
library(tidygraph)

T3K_path <- system.file("extdata", '3k_pan_T_filtered_gene_bc_matrices/GRCh38', package = 'ReSET', mustWork = T)
T3K_data <- ReSET::preprocess_data(T3K_path)
phenograph_results <- ReSET::PhenoGraphVarRes(t(T3K_data@scale.data[T3K_data@var.genes,]), c(100, 50, 40, 30, 20)) 
clusters.aggregate <- cbind(1, phenograph_results$membership)
colnames(clusters.aggregate) = paste('K',1:6)
clustree_res <- ReSET::clustertree_construction(clusters.aggregate, prefix = "K")
print(clustree_res)


T3K_tree <- ReSET::phenograph_construction(clusters.aggregate)
ReSET::plot_graph(T3K_tree$graph)
ReSET::plot_graph(T3K_tree$trimmed)

T4K_path <- system.file("extdata", '4k_pan_T_filtered_gene_bc_matrices/GRCh38', package = 'ReSET', mustWork = T)
T4K_data <- ReSET::preprocess_data(T4K_path)
phenograph_results <- ReSET::PhenoGraphVarRes(t(T4K_data@scale.data[T4K_data@var.genes,]), c(100, 50, 40, 30, 20))
clusters.aggregate <- cbind(1, phenograph_results$membership)
colnames(clusters.aggregate) = paste('K',1:6)
clustree_res <- clustree::clustree(clusters.aggregate, prefix = "K")
print(clustree_res)
T4K_tree <- ReSET::clustertree_construction(clusters.aggregate)
plot_graph(T4K_tree$graph)
plot_graph(T4K_tree$trimmed)
