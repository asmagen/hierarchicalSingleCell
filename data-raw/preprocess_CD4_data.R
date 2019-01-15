# reproducing CD44 replicate 1 dataset (Ciucci et al 2018)
library(Seurat)
library(dplyr)
library(ggplot2)

CD44_1 <- Seurat::Read10X('/Users/szmamie/Documents/repos/robustSingleCell/inst/extdata/CD44Total_1')  # 2211840 x 27998  
CD44_1 <- Seurat::CreateSeuratObject(CD44_1, min.cells = 2, min.genes = 500) # 14313  x 7112
mito_genes <- grepl("^mt-", rownames(CD44_1@data)) 
CD44_1@meta.data$percent_mito <- Matrix::colSums(CD44_1@raw.data[mito_genes, ]) / Matrix::colSums(CD44_1@raw.data)
CD44_1 <- Seurat::FilterCells(CD44_1, subset.names = 'percent_mito', high.thresholds = 0.1) # 14313 x 7021  (reference 7006)
CD44_1 <- Seurat::NormalizeData(CD44_1)
CD44_1 <- Seurat::ScaleData(CD44_1)
CD44_1 <- Seurat::FindVariableGenes(CD44_1, x.low.cutoff = 0.0125, x.high.cutoff = 6, y.cutoff = 0.5)
CD44_1 <- Seurat::RunPCA(CD44_1, pc.genes = CD44_1@var.genes, pcs.compute = 20)
CD44_1 <- Seurat::FindClusters(CD44_1, reduction.type = "pca", dims.use = 1:20, resolution = 2, print.output = 0, save.SNN = F, force.recalc = T)
CD44_1 <- Seurat::SetAllIdent(CD44_1,id = "res.2")
prop.table(table(CD44_1@ident))
CD44_1 <- Seurat::RunTSNE(CD44_1, dims.use = 1:20, do.fast = TRUE)
Seurat::FeaturePlot(CD44_1, c("Cd4","Cd8b1"))
Seurat::TSNEPlot(CD44_1, group.by="res.2")

# subset for CD4 cells
CD44_1 <- Seurat::AddModuleScore(CD44_1, genes.list = list(c('Cd4'), c('Cd8a'), c('Cd8b1')))
cluster_median <- CD44_1@meta.data %>%
  select(res.2, Cluster1, Cluster2, Cluster3) %>%
  rename(cluster = res.2, CD4 = Cluster1, CD8A = Cluster2, CD8b1 = Cluster3) %>%
  group_by(cluster) %>%
  summarize_all(median)
ggplot(data = CD44_1@meta.data) +
  geom_point(aes(x = Cluster1, y = Cluster2, color = res.2), alpha = 0.5, size = 1) +
  theme_classic() +
  xlab('CD4') +
  ylab('CD8b1') +
  geom_text(data = cluster_median, aes(x = CD4, y = CD8A, label = cluster))
ggplot(data = CD44_1@meta.data) +
  geom_point(aes(x = Cluster1, y = Cluster3, color = res.2), alpha = 0.5, size = 1) +
  theme_classic() +
  xlab('CD4') +
  ylab('CD8a') +
  geom_text(data = cluster_median, aes(x = CD4, y = CD8b1, label = cluster))

CD44_CD4_1 <- Seurat::SubsetData(CD44_1, cells.use = which(CD44_1@meta.data$res.2 %in% c(2, 5, 10, 12, 3, 7))) # 2793 (reference was 2782)

# perform highly variable genes on selected cells and cluster
CD44_CD4_1 <- Seurat::FindVariableGenes(CD44_CD4_1, x.low.cutoff = 0.0125, x.high.cutoff = 6, y.cutoff = 0.5)
CD44_CD4_1 <- Seurat::RunPCA(CD44_CD4_1, pc.genes = CD44_CD4_1@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)
CD44_CD4_1 <- Seurat::FindClusters(CD44_CD4_1, reduction.type = "pca", dims.use = 1:20, resolution = 1, print.output = 0, save.SNN = F, force.recalc = T)
CD44_CD4_1 <- Seurat::SetAllIdent(CD44_CD4_1, id = "res.1")
prop.table(table(CD44_CD4_1@ident))

CD44_CD4_1 <- Seurat::RunTSNE(CD44_CD4_1, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(CD44_CD4_1, group.by="res.1")

