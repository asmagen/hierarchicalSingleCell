# reproducing CD44 replicate 1 dataset (Ciucci et al 2018)

library(Seurat)
library(dplyr)
library(ggplot2)
library(ReSET)

CD44_1 <- Seurat::Read10X("~/Desktop/Inbox/CD44Total_1/")  # 2211840 x 27998
CD44_2 <- Seurat::Read10X("~/Desktop/Inbox/CD44Total_2/")  # 737280 x 27998

CD44_1 <- Seurat::CreateSeuratObject(CD44_1, min.cells = 2, min.genes = 500) # 14313  x 7112
CD44_2 <- Seurat::CreateSeuratObject(CD44_2, min.cells = 2, min.genes = 500) # 12845 x 2224

mito_genes <- grepl("^mt-", rownames(CD44_1@data))
CD44_1@meta.data$percent_mito <- Matrix::colSums(CD44_1@raw.data[mito_genes, ]) / Matrix::colSums(CD44_1@raw.data)
CD44_1 <- Seurat::FilterCells(CD44_1, subset.names = 'percent_mito', high.thresholds = 0.1) # 14313 x 7021  (reference 7006)

mito_genes <- grepl("^mt-", rownames(CD44_2@data))
CD44_2@meta.data$percent_mito <- Matrix::colSums(CD44_2@raw.data[mito_genes, ]) / Matrix::colSums(CD44_2@raw.data)
CD44_2 <- Seurat::FilterCells(CD44_2, subset.names = 'percent_mito', high.thresholds = 0.1) # 12845 x 2193  (reference ?)

# cluster CD44_1
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

# subset for CD4 cells for CD44_1
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

# cluster CD44_2
CD44_2 <- Seurat::NormalizeData(CD44_2)

CD44_2 <- Seurat::FindVariableGenes(CD44_2, x.low.cutoff = 0.0125, x.high.cutoff = 6, y.cutoff = 0.5)
CD44_2 <- Seurat::ScaleData(CD44_2)
CD44_2 <- Seurat::RunPCA(CD44_2, pc.genes = CD44_2@var.genes, pcs.compute = 20)
CD44_2 <- Seurat::FindClusters(CD44_2, reduction.type = "pca", dims.use = 1:20, resolution = 2, print.output = 0, save.SNN = F, force.recalc = T) # 13 clusters
CD44_2 <- Seurat::SetAllIdent(CD44_2,id = "res.2")
prop.table(table(CD44_2@ident))
CD44_2 <- Seurat::RunTSNE(CD44_2, dims.use = 1:20, do.fast = TRUE)
Seurat::FeaturePlot(CD44_2, c("Cd4","Cd8b1"))
Seurat::TSNEPlot(CD44_2, group.by="res.2", do.label = TRUE)

# subset for CD4 cells for CD44_2
CD44_2 <- Seurat::AddModuleScore(CD44_2, genes.list = list(c('Cd4'), c('Cd8a'), c('Cd8b1')))
cluster_median <- CD44_2@meta.data %>%
  select(res.2, Cluster1, Cluster2, Cluster3) %>%
  rename(cluster = res.2, CD4 = Cluster1, CD8A = Cluster2, CD8b1 = Cluster3) %>%
  group_by(cluster) %>%
  summarize_all(median)
ggplot(data = CD44_2@meta.data) +
  geom_point(aes(x = Cluster1, y = Cluster2, color = res.2), alpha = 0.5, size = 1) +
  theme_classic() +
  xlab('CD4') +
  ylab('CD8b1') +
  geom_text(data = cluster_median, aes(x = CD4, y = CD8A, label = cluster))
ggplot(data = CD44_2@meta.data) +
  geom_point(aes(x = Cluster1, y = Cluster3, color = res.2), alpha = 0.5, size = 1) +
  theme_classic() +
  xlab('CD4') +
  ylab('CD8a') +
  geom_text(data = cluster_median, aes(x = CD4, y = CD8b1, label = cluster))

# what looks like CD4 cells to me (Billy)
CD44_CD4_2 <- Seurat::SubsetData(CD44_2, cells.use = which(CD44_2@meta.data$res.2 %in% c(1,3,4,7,9,12))) # 963 (reference was 956)

# what Thomas says are CD4 clusters
CD44_CD4_2 <- Seurat::SubsetData(CD44_2, cells.use = which(CD44_2@meta.data$res.2 %in% c(0,2,4,8,9))) # 963 (reference was 956)


# perform highly variable genes on selected cells and cluster CD44_CD4_1
CD44_CD4_1 <- Seurat::FindVariableGenes(CD44_CD4_1, x.low.cutoff = 0.0125, x.high.cutoff = 6, y.cutoff = 0.5)
CD44_CD4_1 <- Seurat::RunPCA(CD44_CD4_1, pc.genes = CD44_CD4_1@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)
CD44_CD4_1 <- Seurat::FindClusters(CD44_CD4_1, reduction.type = "pca", dims.use = 1:20, resolution = 1, print.output = 0, save.SNN = F, force.recalc = T)
CD44_CD4_1 <- Seurat::SetAllIdent(CD44_CD4_1, id = "res.1")
prop.table(table(CD44_CD4_1@ident))

CD44_CD4_1 <- Seurat::RunTSNE(CD44_CD4_1, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(CD44_CD4_1, group.by="res.1", do.label=TRUE)

# perform highly variable genes on selected cells and cluster CD44_CD4_2
CD44_CD4_2 <- Seurat::FindVariableGenes(CD44_CD4_2, x.low.cutoff = 0.0125, x.high.cutoff = 6, y.cutoff = 0.5)
CD44_CD4_2 <- Seurat::RunPCA(CD44_CD4_2, pc.genes = CD44_CD4_2@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)
CD44_CD4_2 <- Seurat::FindClusters(CD44_CD4_2, reduction.type = "pca", dims.use = 1:20, resolution = 1, print.output = 0, save.SNN = F, force.recalc = T)
CD44_CD4_2 <- Seurat::SetAllIdent(CD44_CD4_2, id = "res.1")
prop.table(table(CD44_CD4_2@ident))

CD44_CD4_2 <- Seurat::RunTSNE(CD44_CD4_2, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(CD44_CD4_2, group.by="res.1", do.label=TRUE)




# Find gene signatures for CD44_CD4_1
CD44_CD4_1.mks<-FindAllMarkers(CD44_CD4_1, only.pos = TRUE)
CD44_CD4_1.mks<-CD44_CD4_1.mks[which(CD44_CD4_1.mks$p_val_adj < 0.05), ]
#write.table(CD44_CD4_1.mks, "/data/HBCC_analysis/analysis/kim2/Hackathon/CD44_replicates_alignment/CD44_CD4_1.mks.txt", sep="\t")

# Find gene signatures for CD44_CD4_2
CD44_CD4_2.mks<-FindAllMarkers(CD44_CD4_2, only.pos = TRUE)
CD44_CD4_2.mks<-CD44_CD4_2.mks[which(CD44_CD4_2.mks$p_val_adj < 0.05), ]
write.table(CD44_CD4_2.mks, "/data/HBCC_analysis/analysis/kim2/Hackathon/CD44_replicates_alignment/CD44_CD4_2.mks.txt", sep="\t")

# Overlap of gene signatures
olap.df<-data.frame(matrix(nrow=length(unique(CD44_CD4_1@ident))* length(unique(CD44_CD4_2@ident)), ncol=3))
colnames(olap.df)<-c("CD44_CD4_1", "CD44_CD4_2", "pct.gene.sig.ident")

i=1
for (cl_1 in levels(CD44_CD4_1@ident)){
  print(paste("Replicate 1: Cluster ", cl_1, sep=""))
  rep1_gene_sig<-CD44_CD4_1.mks[which(CD44_CD4_1.mks$cluster == cl_1),]$gene
  rep1_length<-length(rep1_gene_sig)

  for (cl_2 in levels(CD44_CD4_2@ident)){
    print(paste("Replicate 2: Cluster ", cl_2, sep=""))
    rep2_gene_sig<-CD44_CD4_2.mks[which(CD44_CD4_2.mks$cluster == cl_2),]$gene
    rep2_length<-length(rep2_gene_sig)

    overlap<-intersect(rep1_gene_sig, rep2_gene_sig)
    num.overlap<-length(overlap)
    all.genes<-length(unique(c(rep1_gene_sig, rep2_gene_sig)))
    pct.overlap<-100*(num.overlap/all.genes)


    olap.df[i,1]<-cl_1
    olap.df[i,2]<-cl_2
    olap.df[i,3]<-pct.overlap
    i=i+1
  }
}

ggplot(olap.df, aes(x=CD44_CD4_1, y=CD44_CD4_2, fill=pct.gene.sig.ident)) +
  geom_tile() +
  scale_fill_gradientn(colors=rev(c(colorRampPalette(c("red", "white"))(5)[-5],colorRampPalette(c("white", "blue"))(6))))

# CD44_CD4_1: Plot genes in Supp Fig 1C
CD44_CD4_1.avg<-AverageExpression(CD44_CD4_1)
genes<-c("Prdm1","Tbx21", "Id2", "Cxcr6", "Bcl6", "Pdcd1", "Cxcr5", "Tcf7", "Id3", "Bcl2", "Ccr7", "Il7r", "Foxp3")
CD44_CD4_1.avg2<-CD44_CD4_1.avg[genes,]
CD44_CD4_1.scl<-t(scale(t(CD44_CD4_1.avg2)))
library(pheatmap)
cols<-rev(RColorBrewer::brewer.pal(9, "RdBu"))
pheatmap(CD44_CD4_1.scl, color = cols, cluster_rows = FALSE)

# CD44_CD4_2: Plot genes in Supp Fig 1C
CD44_CD4_2.avg<-AverageExpression(CD44_CD4_2)
CD44_CD4_2.avg2<-CD44_CD4_2.avg[genes,]
CD44_CD4_2.scl<-t(scale(t(CD44_CD4_2.avg2)))
pheatmap(CD44_CD4_2.scl, color = cols, cluster_rows = FALSE)

# Rename cells based on expression of genes in Supp Fig 1C&D
df.1<-data.frame(cl.mem=c(4,2,5,3,1,0),
                 cl.name=c("Th1","Tfh","Tcmp","Tmem","Treg","Rep1_Unk1"),
                 stringsAsFactors = FALSE)
df.2<-data.frame(cl.mem=c(4,2,3,1,0,5),
                 cl.name=c("Th1","Tfh","Tcmp","Treg","Rep2_Unk1","Rep2_Unk2"),
                 stringsAsFactors = FALSE)

CD44_CD4_1@meta.data$ID<-plyr::mapvalues(CD44_CD4_1@meta.data$res.1,
                                         from = df.1$cl.mem,
                                         to = df.1$cl.name)
CD44_CD4_1<-SetAllIdent(CD44_CD4_1, id = "ID")
CD44_CD4_2@meta.data$ID<-plyr::mapvalues(CD44_CD4_2@meta.data$res.1,
                                         from = df.2$cl.mem,
                                         to = df.2$cl.name)
CD44_CD4_2<-SetAllIdent(CD44_CD4_2, id = "ID")


# recalculate average expression
CD44_CD4_1.avg<-AverageExpression(CD44_CD4_1)
CD44_CD4_2.avg<-AverageExpression(CD44_CD4_2)

CD44_CD4_1.avg2<-CD44_CD4_1.avg[genes,]
CD44_CD4_1.scl<-t(scale(t(CD44_CD4_1.avg2)))
pheatmap(CD44_CD4_1.scl, color = cols, cluster_rows = FALSE)

CD44_CD4_2.avg2<-CD44_CD4_2.avg[genes,]
CD44_CD4_2.scl<-t(scale(t(CD44_CD4_2.avg2)))
pheatmap(CD44_CD4_2.scl, color = cols, cluster_rows = FALSE)


# hierarchy construction
CD44_CD4_1_hierarchy <- ReSET::construct_hierarchy(t(CD44_CD4_1@scale.data[CD44_CD4_1@var.genes, ]),
                                              membership = CD44_CD4_1@ident)
stats:::plot.dendrogram(CD44_CD4_1_hierarchy$dend)
CD44_CD4_1_tree <- with(CD44_CD4_1_hierarchy, as_binary_tree(hclust, data, membership))

CD44_CD4_2_hierarchy <- ReSET::construct_hierarchy(t(CD44_CD4_2@scale.data[CD44_CD4_2@var.genes, ]),
                                                   membership = CD44_CD4_2@ident)
stats:::plot.dendrogram(CD44_CD4_2_hierarchy$dend)
CD44_CD4_2_tree <- with(CD44_CD4_2_hierarchy, as_binary_tree(hclust, data, membership))

par(mfrow=c(1,2))
plot(CD44_CD4_1_tree$tree)
plot(CD44_CD4_2_tree$tree)

cost_matrix <- ReSET::cor_cost_matrix(CD44_CD4_1_tree$summary_stats, CD44_CD4_2_tree$summary_stats)

align_obj <- ReSET::align(CD44_CD4_1_tree$tree, CD44_CD4_2_tree$tree, cost_matrix)
plot(align_obj$tree)
