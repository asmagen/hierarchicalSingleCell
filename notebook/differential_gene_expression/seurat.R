library(Seurat)
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite('BiocManager')
BiocManager::install("biomaRt", suppressUpdates = TRUE)
library(biomaRt)

#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
#install.packages('cellrangerRkit')
#library(cellrangerRkit)

#gbm <- load_cellranger_matrix(file.path(environment$data.path,dataset), genome=genome)
#measurements = as.matrix(exprs(gbm)[,colSums(exprs(gbm)>0)>=min.genes.per.cell])
#rownames(measurements) = fData(gbm)$symbol


t4 <- as.data.frame(fread('t_4k_matrix.csv'))
t3 <- as.data.frame(fread('t_3k_matrix.csv'))

genes <- t4$V1
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
G_list <- as.data.table(G_list)
genes <- as.data.table(genes)
genes <- merge(genes, G_list, all=T, by.x="genes", by.y="ensembl_gene_id")
genes <- na.omit(genes)
t4 <- merge(t4, genes, by.x='V1', by.y='genes', all=T)
t4 <- na.omit(t4)
t4 <- as.data.table(t4)
t4 <- t4[,c(4540, 2:4539)]
t4 <- t4[!duplicated(t4[,c(external_gene_name)] ),]
fwrite(t4, 't_4k_genesymbol.csv', sep='\t', col.names = T, row.names = F, quote=F)
t4 <- as.data.frame(fread('t_4k_genesymbol.csv'))
t4 <- as.data.frame(t4)
row.names(t4) <- t4$external_gene_name
t4 <- t4[1:33201, 2:4539]
t4 <- CreateSeuratObject(raw.data = t4, project = "t4")
t4@meta.data$sample <- "t4"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = t4@data), value = TRUE)
percent.mito <- Matrix::colSums(t4@raw.data[mito.genes, ])/Matrix::colSums(t4@raw.data)
rpl.genes <- grep(pattern = "^RPL", x = rownames(x = t4@data), value = TRUE)
rps.genes <- grep(pattern = "^RPS", x = rownames(x = t4@data), value = TRUE)
percent.rps <- Matrix::colSums(t4@raw.data[rps.genes, ])/Matrix::colSums(t4@raw.data)
percent.rpl <- Matrix::colSums(t4@raw.data[rpl.genes, ])/Matrix::colSums(t4@raw.data)

t4 <- AddMetaData(object = t4, metadata = percent.mito, col.name = "percent.mito")
t4 <- AddMetaData(object = t4, metadata = percent.rpl, col.name = "percent.rpl")
t4 <- AddMetaData(object = t4, metadata = percent.rps, col.name = "percent.rps")

#VlnPlot(object = t4, features.plot = c("nGene", "nUMI", "percent.mito", "percent.rpl", "percent.rps"), nCol = 5)

#GenePlot(object = t4, gene1 = "nUMI", gene2 = "percent.mito" )
#GenePlot(object = t4, gene1 = "nUMI", gene2 = "nGene")

# percent mito threshold: 0 to 0.2
# ngenes : 500 to 3000

t4 <- FilterCells(object = t4, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, 0), high.thresholds = c(3000, 0.2))
t4 <- NormalizeData(t4)
t4 <- ScaleData(t4, display.progress = F, vars.to.regress = c('nGene', 'percent.mito', 'percent.rpl', 'percent.rps'))
t4 <- FindVariableGenes(object = t4, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
t4 <- FindClusters(object = t4, reduction.type = "pca", dims.use = 1:10, resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
t4 <- FindClusters(object = t4, reduction.type = "pca", dims.use = 1:10, resolution = 0.1, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
t4 <- FindClusters(object = t4, reduction.type = "pca", dims.use = 1:10, resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
t4 <- FindClusters(object = t4, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)

t4 <- SetAllIdent(t4, id='res.0.6')
#t4.0 <- FindMarkers(t4, ident.1 = '0')
#View(t4.0.vs.5)

t4 <- SetAllIdent(t4, id="res.0.6")
t4 <- RunPCA(object = t4, pcs.print = 1:5, genes.print = 10)
#PCAPlot(t4)
t4 <- RunTSNE(object = t4, dims.use = 1:8, do.fast = TRUE)

Th1 <- c("CD4", "CXCR3", "CCR5", "IL12A")
Th2 <- c("CD4", "CCR4", "IL4", "IL5", "GATA3")
Th9 <- c("CD4", "CCR3", "CCR6", "IL9", "SFPI1")
Th17 <- c("CD4", "CCR6", "CCR4", "KLRB1", "IL17")
Th22 <- c("CD4", "CCR10", "CCR4", "CCR6", "IL22")
Tfh <- c("CD4", "CXCR5", "CD40L", "ICOS", "IL21", "BCL6")

cytotoxic.cd8.t <- c("CD8", "IFNG", "PRF1", "GZMA")


t4 <- SetAllIdent(t4, id='res.0.6')
VlnPlot(t4, c(Th1))
FeaturePlot(t4, c("CD4", "CXCR3", "CCR5", "IL12A", "IFNG", "TBX21"), cols.use = c('green', 'red'))
VlnPlot(t4, c(Th1, Th2, Th9, Th17, Th22, Tfh))

# t3----
t3 <- as.data.frame(fread('t_3k_matrix.csv'))
genes <- t3$V1
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
G_list <- as.data.table(G_list)
genes <- as.data.table(genes)
genes <- merge(genes, G_list, all=T, by.x="genes", by.y="ensembl_gene_id")
genes <- na.omit(genes)
t3 <- merge(t3, genes, by.x='V1', by.y='genes', all=T)
t3 <- na.omit(t3)
t3 <- as.data.table(t3)
t3 <- t3[!duplicated(t3[,c(external_gene_name)] ),]
t3 <- t3[,c(3557, 2:3556)]
fwrite(t3, 't_3k_genesymbol.csv', sep='\t', col.names = T, row.names = F, quote=F)
t3 <- as.data.frame(t3)
row.names(t3) <- t3$external_gene_name
t3 <- t3[1:33202, 2:3556]
t3 <- CreateSeuratObject(raw.data = t3, project = "t3")
t3@meta.data$sample <- "t3"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = t3@data), value = TRUE)
percent.mito <- Matrix::colSums(t3@raw.data[mito.genes, ])/Matrix::colSums(t3@raw.data)
rpl.genes <- grep(pattern = "^RPL", x = rownames(x = t3@data), value = TRUE)
rps.genes <- grep(pattern = "^RPS", x = rownames(x = t3@data), value = TRUE)
percent.rps <- Matrix::colSums(t3@raw.data[rps.genes, ])/Matrix::colSums(t3@raw.data)
percent.rpl <- Matrix::colSums(t3@raw.data[rpl.genes, ])/Matrix::colSums(t3@raw.data)
t3 <- AddMetaData(object = t3, metadata = percent.mito, col.name = "percent.mito")
t3 <- AddMetaData(object = t3, metadata = percent.rpl, col.name = "percent.rpl")
t3 <- AddMetaData(object = t3, metadata = percent.rps, col.name = "percent.rps")
t3 <- FilterCells(object = t3, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, 0), high.thresholds = c(3000, 0.2))
t3 <- NormalizeData(t3)
t3 <- ScaleData(t3, display.progress = F)
t3 <- RunPCA(object = t3, pcs.print = 1:5, genes.print = 5)
t3 <- FindVariableGenes(object = t3, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
t3 <- FindClusters(object = t3, reduction.type = "pca", dims.use = 1:10, resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
t3 <- FindClusters(object = t3, reduction.type = "pca", dims.use = 1:10, resolution = 0.1, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
t3 <- FindClusters(object = t3, reduction.type = "pca", dims.use = 1:10, resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
t3 <- FindClusters(object = t3, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
t3 <- SetAllIdent(t3, id='res.0.6')
t3 <- RunPCA(object = t3, pcs.print = 1:5, genes.print = 5)
t3 <- RunTSNE(object = t3, dims.use = 1:8, do.fast = TRUE)

save(t3, file='t3.Seurat.RData')