library(Seurat)
library(dplyr)
library(reshape2)

# Load the PBMC dataset
k3.data <- Read10X(data.dir = "/ldata/single/Hackathon/3k_filtered_gene_bc_matrices/GRCh38/")
k4.data <- Read10X(data.dir = "/ldata/single/Hackathon/4k_filtered_gene_bc_matrices/GRCh38/")
# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size
sparse.size <- object.size(x = pbmc.data)
sparse.size
dense.size / sparse.size

pbmc3 <- CreateSeuratObject(raw.data = k3.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")
pbmc4 <- CreateSeuratObject(raw.data = k4.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc3@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc3@raw.data[mito.genes, ]) / Matrix::colSums(pbmc3@raw.data)
# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
pbmc3 <- AddMetaData(object = pbmc3, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc3, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc4@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc4@raw.data[mito.genes, ]) / Matrix::colSums(pbmc4@raw.data)

pbmc4 <- AddMetaData(object = pbmc4, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc4, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


# GenePlot is typically used to visualize gene-gene relationships, but can be used for anything 
# calculated by the object, i.e. columns in object@meta.data, PC scores etc.
# Since there is a rare subset of cells with an outlier level of high mitochondrial percentage
# and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than 200
# Note that low.thresholds and high.thresholds are used to define a 'gate'.
# -Inf and Inf should be used if you don't want a lower or upper threshold.
pbmc3 <- FilterCells(object = pbmc3, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
pbmc4 <- FilterCells(object = pbmc4, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

pbmc3 <- NormalizeData(object = pbmc3, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc4 <- NormalizeData(object = pbmc4, normalization.method = "LogNormalize", scale.factor = 1e4)

### Detection of variable genes across the single cells

pbmc3 <- FindVariableGenes(object = pbmc3, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc4 <- FindVariableGenes(object = pbmc4, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

pbmc3 <- ScaleData(object = pbmc3, vars.to.regress = c("nUMI", "percent.mito"))
pbmc4 <- ScaleData(object = pbmc4, vars.to.regress = c("nUMI", "percent.mito"))

pbmc3 <- RunPCA(object = pbmc3, pc.genes = pbmc3@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc4 <- RunPCA(object = pbmc4, pc.genes = pbmc4@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)


PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)

JackStrawPlot(object = pbmc, PCs = 1:12)

PCElbowPlot(object = pbmc3)
PCElbowPlot(object = pbmc4)

pbmc3 <- FindClusters(object = pbmc3, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
pbmc4 <- FindClusters(object = pbmc4, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)

pbmc3 <- RunTSNE(object = pbmc3, dims.use = 1:10)
pbmc4 <- RunTSNE(object = pbmc4, dims.use = 1:10)

TSNEPlot(object = pbmc3)
TSNEPlot(object = pbmc4)

saveRDS(pbmc3, file="/ldata/single/Hackathon/pbmc3.rds")
saveRDS(pbmc4, file="/ldata/single/Hackathon/pbmc4.rds")
#saveRDS(pbmc, file = "~/Projects/datasets/pbmc3k/pbmc_tutorial.rds")

# find all markers of cluster 1
#cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
#print(x = head(x = cluster1.markers, n = 5))
# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
#print(x = head(x = cluster5.markers, n = 5))
# find markers for every cluster compared to all remaining cells, report only the positive ones
#pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

###############################################################################################################
# Find intersecting variable genes
gs<-intersect(pbmc3@var.genes,pbmc3@var.genes)
length(gs)

############################## Manual average expression calculation ######################################
avg_exp<-function(seur_obj){
  h<-as.data.frame(as.matrix(seur_obj@data))
  h$gene<-rownames(h)
  h<-melt(h,id="gene")
  #
  t<-seur_obj@meta.data
  t2<-cbind(t,rownames(t))
  head(t2)
  #
  colnames(t2)<-c("nGene","nUMI","orig.ident","percent.mito","res.0.6","cell")
  h2<-h %>% mutate(cell=variable) %>% left_join(t2[,c("res.0.6","cell")],by="cell")
  head(h2)

  h3<-h2 %>% group_by(gene, res.0.6) %>% dplyr::summarise(meanExp=mean(value)) %>% ungroup
  h3<-as.data.frame(h3)
  #
  h4<-dcast(h3, gene ~ res.0.6, value.var="meanExp")
  head(h4)
  rownames(h4)<-h4$gene
  h4<-as.matrix(h4[,-1])

  return(h4)
}
p3<-avg_exp(pbmc3)
p4<-avg_exp(pbmc4)
############################## Manual average z-score calculation ############################
avg_z<-function(seur_obj){
  h<-as.data.frame(as.matrix(seur_obj@scale.data))
  h$gene<-rownames(h)
  h<-melt(h,id="gene")
  #
  t<-seur_obj@meta.data
  t2<-cbind(t,rownames(t))
  head(t2)
  #
  colnames(t2)<-c("nGene","nUMI","orig.ident","percent.mito","res.0.6","cell")
  h2<-h %>% mutate(cell=variable) %>% left_join(t2[,c("res.0.6","cell")],by="cell")
  head(h2)
  
  h3<-h2 %>% group_by(gene, res.0.6) %>% dplyr::summarise(meanExp=mean(value)) %>% ungroup
  h3<-as.data.frame(h3)
  #
  h4<-dcast(h3, gene ~ res.0.6, value.var="meanExp")
  head(h4)
  rownames(h4)<-h4$gene
  h4<-as.matrix(h4[,-1])
  
  return(h4)
}
z3<-avg_z(pbmc3)
z4<-avg_z(pbmc4)

######################################### heirarchical clustering with hclust ################
h_clust<-function(avg_exp_output){
  avg_exp_output<-t(scale(t(avg_exp_output),scale=T))
  dim(avg_exp_output)
  #
  h.mak<-avg_exp_output[pbmc4@var.genes,]
  dim(h.mak)
  d<-dist(as.matrix(t(h.mak)))
  hcc<-hclust(d,method="ward.D")
  hcc<-as.dendrogram(hcc)
  return(hcc)
}
p3_hc<-h_clust(p3)
p4_hc<-h_clust(p4)
plot(p3_hc)
plot(p4_hc)
############################# Average Expression and Pearson Correlation Heatmap ######################
p3<-AverageExpression(object=pbmc3)
p4<-AverageExpression(object=pbmc4)
get_cor_df<-function(emat1, emat2){
  genelist<-intersect(pbmc3@var.genes,pbmc4@var.genes) #only consider co-expressed markers
  print(length(genelist))
  df1<-emat1[genelist,]
  df2<-emat2[genelist,]
  
  #df1<-log1p(df1)
  #df2<-log1p(df2)
  #try this outside
  num<-as.numeric(dim(df1)[2]*dim(df2)[2])
  df.cor<-data.frame(matrix(ncol=3,nrow=num))
  
  k=0
  for(i in colnames(df2)){
    for(j in colnames(df1)){
      k=k+1
      df.cor[k,]<-c(i,j,cor(df2[,i],df1[,j],use="complete.obs"))
      #print(c(i,j,cor(df2[[i]],df1[[j]])))
      #print(k)
    }
  }
  df.cor$X3<-as.numeric(df.cor$X3)
  return(df.cor)
}

comp.n<-get_cor_df(p3,p4)
ggplot(data=comp.n,aes(x=X2,y=X1,fill=X3))+geom_tile() + 
  scale_fill_gradientn(colors=rev(c(colorRampPalette(c("royalblue4", "white"))(5)[-5],colorRampPalette(c("white", "gold2"))(6)))) +
  theme(axis.text.x=element_text(angle = 45,hjust=1,vjust=1), axis.ticks.x=element_blank(),axis.title.x=element_blank())+ xlab("")+
  ylab("") + theme(legend.position="bottom",legend.title=element_blank()) + ggtitle("")

comp.n2<-get_cor_df(z3,z4)
ggplot(data=comp.n2,aes(x=X2,y=X1,fill=X3))+geom_tile() + 
  scale_fill_gradientn(colors=rev(c(colorRampPalette(c("royalblue4", "white"))(5)[-5],colorRampPalette(c("white", "gold2"))(6)))) +
  theme(axis.text.x=element_text(angle = 45,hjust=1,vjust=1), axis.ticks.x=element_blank(),axis.title.x=element_blank())+ xlab("")+
  ylab("") + theme(legend.position="bottom",legend.title=element_blank()) + ggtitle("")
############################################# DE ############################################

output<-function(input){
  function(){
    
    
  }
  
  
}
############################################# Jaccard ########################################
# input: 2 matrices of n (diff. expressed genes) x 2 --> 1 for each experiment
inputA<-data.frame(gene=c("a","b","c","d","e","f","g"),node=c(121,121,121,121,122,122,122))
inputB<-data.frame(gene=c("a","b","c","d","e","f","g"),node=c(123,123,124,124,124,124,124))
get_jaccard<-function(inputA, inputB){
  splitA<-split(inputA,inputA$node)
  splitB<-split(inputB,inputB$node)
  num<-as.numeric(length(splitA)*length(splitB))
  jc<-data.frame(matrix(ncol=3,nrow=num))
  k=0
  for (i in seq(length(splitA))){
    v1<-splitA[[i]]$gene
    for (j in seq(length(splitB))){
      k=k+1
      print
      v2<-splitB[[j]]$gene
      jaccard_calc<-length(intersect(v1,v2))/length(union(v1,v2))
      jc[k,]<-c(names(splitA)[i],names(splitB)[j],jaccard_calc)
    }
  }
  jc$X3<-as.numeric(jc$X3)
  jc2<-matrix(jc$X3,nrow=length(unique(jc$X1)),ncol=length(unique(jc$X2)), byrow = TRUE)
  colnames(jc)<-c("node_exp1","node_exp2","j_indx")
  rownames(jc2)<-names(splitA)
  colnames(jc2)<-names(splitB)
  return(list(jc=jc,jc2=jc2))
}
j<-get_jaccard(inputA,inputB)

library(gplots)
heatmap.2(j$jc2, col=bluered, trace="none")


##############################################################################################
ggplot(data=pbmc3@meta.data,aes(x=res.0.6,y=nGene))+geom_violin()+
  scale_y_continuous(limits=c(0,2500))
ggplot(data=pbmc4@meta.data,aes(x=res.0.6,y=nGene))+geom_violin()+
  scale_y_continuous(limits=c(0,2500))

###########################################################################################
ggplot(data=pbmc3@meta.data,aes(x=res.0.6,y=nUMI))+geom_violin()+
  scale_y_continuous(limits=c(0,10000))
ggplot(data=pbmc4@meta.data,aes(x=res.0.6,y=nUMI))+geom_violin()+
  scale_y_continuous(limits=c(0,10000))
###########################################################################################