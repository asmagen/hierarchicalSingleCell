library(Seurat)
library(dplyr)

# Load the PBMC dataset, filter and normalize
pbmc.data <- Read10X(data.dir = "/Users/rhodesct/Documents/hackathon/4k Pan T Cells from a Healthy Donor/filtered_gene_bc_matrices/GRCh38/")

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ]) / Matrix::colSums(pbmc@raw.data)

pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their correlation 
# with the calculated components. Though we don't use this further here, it can be used to identify markers that 
# are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection. 
# The results of the projected PCA can be explored by setting use.full=T in the functions above
# pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
# 
# PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# 
# PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
# 
# pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
# 
# JackStrawPlot(object = pbmc, PCs = 1:12)
# 
# PCElbowPlot(object = pbmc)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10)

TSNEPlot(object = pbmc)

#find markers for all groups
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top20 = pbmc.markers %>% group_by(cluster) %>% top_n(-20, p_val_adj) %>% select(c(cluster, p_val_adj, gene))


table(pbmc@ident)
prop.table(x = table(pbmc@ident))
# WhichCells(object = pbmc, ident = "0")
cluster.averages <- AverageExpression(object = pbmc, return.seurat = TRUE, show.progress = FALSE)
cluster_mat <- AverageExpression(object = pbmc)

# filter matrix by top 20 genes per cluster
head(cluster_mat)
include_list <- top20$gene
filtered_mat = cluster_mat[include_list, ]


