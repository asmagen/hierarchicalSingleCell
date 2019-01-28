#' Preprocess 10X dataset 
#' @param path The path to the dataset
#' @param min.cells The minimum number of cells for each gene
#' @param min.genes The minimum number of genes for each cell
#' @param nGene_thresh The quantile threshold for multilplets based on nGene
#' @param percent_mito_thresh  The quantile threshold for multilplets based on mitochondrial gene expression percentage
#' @return A Seurat object with data normalized, scaled and PCA performed
#' @export
#' @author Meng Wang
preprocess_data <- function(path, min.cells = 3, min.genes = 200, nGene_thresh = NULL,
                            percent_mito_thresh = 0.98, mito_thresh = NULL,
                            x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1, regress = F) {
  data <- Seurat::Read10X(path)
  data <- Seurat::CreateSeuratObject(data, min.cells = min.cells, min.genes = min.genes)
  mito_genes <- grepl("^MT-", rownames(data@data))
  data@meta.data$percent_mito <- Matrix::colSums(data@raw.data[mito_genes, ]) / Matrix::colSums(data@raw.data)
  if (!is.null(nGene_thresh)) {
    nGene_high_thresh <- quantile(data@meta.data$nGene, nGene_thresh)
    data <- Seurat::FilterCells(data, subset.names = 'nGene', high.thresholds = nGene_high_thresh)
  }
  if (!is.null(mito_thresh)) percent_mito_high_thresh <- mito_thresh
  else percent_mito_high_thresh <- quantile(data@meta.data$percent_mito, percent_mito_thresh)
  
  data <- Seurat::FilterCells(data, subset.names = 'percent_mito', 
                              high.thresholds = percent_mito_high_thresh)
  
  data <- Seurat::NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  data <- Seurat::FindVariableGenes(data) #mean.function = Seurat::ExpMean, dispersion.function = Seurat::LogVMR, 
                                    #x.low.cutoff = x.low.cutoff, 
                                    #x.high.cutoff = x.high.cutoff, 
                                    #y.cutoff = y.cutoff)
  if (regress) {
    data <- Seurat::ScaleData(data, vars.to.regress = c('nUMI', 'percent_mito'))
  } else {
    data <- Seurat::ScaleData(data)
  }
  data <- Seurat::RunPCA(data, pc.genes = data@var.genes, do.print = F)
  return(data)
}

#' Get highly variable genes
#' @param filtered.normalized Filtered and normalized data matrix
#' @param dir_path Path to the directory where the figure is output
#' @param min.mean The minimum mean for filtering criterion
#' @param min.frac.cells The fraction of cells that are not zero for genes
#' @param min.dispersion.scaled The minimum scaled dispersion for the genes
#' @return The names of the genes that are defined as highly variable based on the criterion
#' @export
get_HVG <- function(filtered.normalized, dir_path, min.mean, min.frac.cells, min.dispersion.scaled,
                    remove_RP = TRUE) {
  if(remove_RP) {
    is_RP <- grepl('^RP[0-9]+[-]', rownames(T3K$normalized_data), fixed = F)
    filtered.normalized <- filtered.normalized[!is_RP,]
  }
  
  means <- apply(filtered.normalized, 1, function(v) log(mean(exp(v) - 1) + 1))
  frac.cells <- rowSums(filtered.normalized > 0)
  vars <- apply(filtered.normalized, 1, function(v) log(var(exp(v) - 1) + 1))
  dispersion <- apply(filtered.normalized, 1, function(v) log(var(exp(v) - 1) / mean(exp(v) - 1)))
  dispersion[is.na(x = dispersion)] <- 0
  means[is.na(x = means)] <- 0
  
  
  pdf(file.path(dir_path, 'VariableGenes.pdf'))
  plot.data <- data.frame(gene = names(means), means, dispersion)
  print(ggplot(plot.data, aes(x=means, y=dispersion, label = gene)) + geom_text(check_overlap = TRUE,size=2) + theme_classic())
  smoothScatter(means,dispersion)
  dev.off()
  
  num.bin <- 20
  bins = cut(x = means, breaks = num.bin)
  names(x = bins) = names(x = means)
  mean_y = tapply(dispersion,bins,mean)
  sd_y = tapply(dispersion,bins,sd)
  dispersion.scaled = (dispersion - mean_y[as.numeric(x = bins)]) / sd_y[as.numeric(x = bins)]
  dispersion.scaled[is.na(x = dispersion.scaled)] = 0
  names(x = dispersion.scaled) = names(x = means)
  
  criterion = means >= min.mean & frac.cells >= min.frac.cells & dispersion.scaled >= min.dispersion.scaled
  HVG = names(means)[criterion]
  return(HVG)
}

#' PhenoGraph Clustering with Various Resolution
#'
#' Run RPhenograph with different resolution and save the membership as a data
#' frame with column names as the k parameter.
#'
#' @param data A input data matrix for Rphenograph.
#' @param K A numeric vector of resolutions.
#' @return A list containing a data matrix of membership at each K and a
#' @return A list containing a data matrix of membership at each K and a
#' vector of clustering modularity
#' @export
#' @examples
#' iris_unique = unique(iris) # Remove duplicates
#' data = as.matrix(iris_unique[,1:4])
#' K = seq(10, 50, length.out=5)
#' res = PhenoGraphVarRes(data, K=K)
#' res
PhenoGraphVarRes <- function(data, K) {
  membership.all = c()
  modularity.all = c()
  numClusters.all = c()
  for (k in K) {
    clusters = Rphenograph::Rphenograph(data, k=k)
    clusters.modularity = igraph::modularity(clusters[[2]])
    clusters.membership = igraph::membership(clusters[[2]])
    numClusters = length(unique(clusters.membership))
    modularity.all = c(modularity.all, clusters.modularity)
    numClusters.all = c(numClusters.all, numClusters)
    membership.all = cbind(membership.all, clusters.membership)
  }
  colnames(membership.all) = K
  return(list(k=K, modularity=modularity.all, numCluster=numClusters.all,
              membership=membership.all))
}