#' return the output of cutree on an hclust object generated from the first
#' ten principal components of the data slot of a seurat object.
#' @param seurat a seurat-class objec
#' @param n an integer (default 16) of the maximum number of splits shown
#' @return a matrix with group memberships from 1 to n
makeCutTree <- function(seurat, n=16) {
  exp <- seurat@data
  exp <- exp[rowSums(exp)>0,]
  exp <- exp[,colSums(exp)>0]
  pca = prcomp(t(exp), retx = TRUE, center = T, scale. = T)
  PCA = t(pca$x)[seq(10),]
  d <- dist(t(PCA))
  tree <- hclust(d, method = "ward.D")
  cutree(tree, k=1:n)  
}
