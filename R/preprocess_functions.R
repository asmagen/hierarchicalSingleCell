#' Preprocess 10X dataset 
#' @param path The path to the dataset
#' @param min.cells The minimum number of cells for each gene
#' @param min.genes The minimum number of genes for each cell
#' @param nGene_thresh The quantile threshold for multilplets based on nGene
#' @param percent_mito_thresh  The quantile threshold for multilplets based on mitochondrial gene expression percentage
#' @return A Seurat object with data normalized, scaled and PCA performed
#' @export
#' @author Meng Wang
preprocess_data <- function(path, min.cells = 3, min.genes = 200, nGene_thresh = 0.98,
                            percent_mito_thresh = 0.98) {
  data <- Seurat::Read10X(path)
  data <- Seurat::CreateSeuratObject(data, min.cells = min.cells, min.genes = min.genes)
  mito_genes <- grepl("^MT-", rownames(data@data))
  data@meta.data$percent_mito <- Matrix::colSums(data@raw.data[mito_genes, ]) / Matrix::colSums(data@raw.data)
  nGene_high_thresh <- quantile(data@meta.data$nGene, nGene_thresh)
  percent_mito_high_thresh <- quantile(data@meta.data$percent_mito, percent_mito_thresh)
  data <- Seurat::FilterCells(data, subset.names = c('nGene', 'percent_mito'), 
                              high.thresholds = c(nGene_high_thresh, percent_mito_high_thresh))
  data <- Seurat::NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  data <- Seurat::FindVariableGenes(data, mean.function = Seurat::ExpMean,
                                    dispersion.function = Seurat::LogVMR, 
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  data <- Seurat::ScaleData(data, vars.to.regress = c('nUMI', 'percent_mito'))
  data <- Seurat::RunPCA(data, pc.genes = data@var.genes, do.print = F)
  return(data)
}
