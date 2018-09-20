#' Preprocess two 10X Genomics pan T cell datasets and one fibroblast
#' reprogamming dataset. Peform QC and normalize.
#' @param input_data the input data.frame with genes in row.names and
#' cell names as columns
#' @return a data.frame with QC filtered cells, genes, and normalized counts

#' To run on the three test datasets, load raw data:
#' library(dyno)
#' library(SingleCellExperiment)
#' library(scater)
#' options(stringsAsFactors = FALSE)
#' library(scran)
#' library(mvoutlier)
#' library(Seurat)
#' data("fibroblast_reprogramming_treutlein")
#' fibro <- as.data.frame(t(fibroblast_reprogramming_treutlein$counts))
#' load('t3.Seurat.RData')
#' load('t4.Seurat.RData')
#' t3.data <- as.data.frame(as.matrix(t3@raw.data))
#' t4.data <- as.data.frame(as.matrix(t4@raw.data))

qc <- function(input_data){
  x <- input_data
  
  # remove ribosomal protein genes
  rp.genes <- row.names(x)[grep('^RP', row.names(x))]
  
  x <- subset(x, !(row.names(x) %in% c(rp.genes)))
  
  # convert to sce
  cellinfo <- data.frame(Cell=c(names(x)), row.names = names(x))
  x <- SingleCellExperiment(assays = list(counts = as.matrix(x)), colData=cellinfo)
  
  # calculate cpm
  exprs(x) <- log2(calculateCPM(x, use_size_factors = FALSE) + 1) 
  
  # calculate qc metrics
  isSpike(x, "MT") <- grepl("^MT-", rownames(x))
  x <- calculateQCMetrics(x,feature_controls = list(MT = isSpike(x, "MT") ))
  
  # remove genes not expressed in any cell
  keep_feature <- rowSums(counts(x) > 0 ) > 0
  x <- x[keep_feature,]
  
  # remove cells with high mitochondrial gene expression
  mito.keep <- !(isOutlier(x$pct_counts_MT, nmads=4, type="higher"))
  
  # remove cells with low library size
  libsize.keep <- !isOutlier(x$total_counts, nmads=4, type="lower", log=TRUE)
  
  # keep cells with sufficient molecules counted and with non-outlier number of reads in MT genes
  x$use <- (libsize.keep & mito.keep)

  x <- computeSumFactors(x)
  x <- normalize(x,exprs_values="counts")
  
  # subset data to usable cells
  x.filt <- as.data.frame(exprs(x))
  x.filt <- x.filt[,x@colData$use]
  
  return(x.filt)
}

#' Perform PCA on normalized counts.
#' @param input_data the normalized counts to run PCA on.
#' @return PCA components.

reset_pca <- function(input_data){
  x <- prcomp(input_data)
  return(x$x)
}
