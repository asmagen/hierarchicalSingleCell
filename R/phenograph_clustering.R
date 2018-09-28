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
PhenoGraphVarRes = function(data, K) {
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


#' Construct Hierarchy based on Cluster Mean
#'
#' Construct cluster hierarchy based on the mean of each cluster.
#'
#' @param data A data matrix to compute cluster mean..
#' @param membership A vector of cluster membership.
#' @return A hclust object.
#' @export
ConstructHierarchy = function(data, membership, func) {
  cluster.mean = c()
  for (cluster in unique(membership)) {
    cluster.data = data[membership == cluster,]
    cluster.mean = rbind(cluster.mean, apply(cluster.data, 2, func))
  }
  rownames(cluster.mean) = as.character(unique(membership))
  dist = stats::dist(cluster.mean, method='euclidean')
  hclust = stats::hclust(dist, method='complete')
  dend = as.dendrogram(hclust)
  dend_k = dendextend::find_k(dend)
  dend = dendextend::color_branches(dend, dend_k$k)
  dend2 = dendextend::seriate_dendrogram(dend, dist)
  ret = list(data=cluster.mean, dend=dend2)
  class(ret) = "hierarchy"
  return(ret)
}

#' Plot Hierarchy Constructed
#'
#' @param x A hierarchy object.
#' @return A heatmap.2 object
#' @export
plot.hierarchy = function(x, ha_row=NULL, show_row_names=T, show_col_names=F) {
  d.marker = stats::dist(t(x$data), method='euclidean')
  hc.marker = stats::hclust(d.marker, method='ward.D2')
  dend = stats::as.dendrogram(hc.marker)
  dend_k = dendextend::find_k(dend)
  dend = dendextend::color_branches(dend, dend_k$k)
  dend2 = dendextend::seriate_dendrogram(dend, d.marker)
  if (!is.null(ha_row)) {
    ComplexHeatmap::Heatmap(x$data, cluster_columns=dend2, cluster_rows=x$dend,
                            col=viridis::viridis(200), show_row_names=show_row_names,
                            show_column_names = show_col_names) + ha_row
  } else {
    ComplexHeatmap::Heatmap(x$data, cluster_columns=dend2, cluster_rows=x$dend,
                            col=viridis::viridis(200), show_row_names = show_row_names,
                            show_column_names = show_col_names)
  }
  
}

#' Test for significant PCs by Permutation Parallel Analysis
#'
#' Estimate a number of significant principal components from a permutation test
#'
#' @param dat A data matrix containing the raw data (each column should be a feature).
#' @param B The number of permutations.
#' @param threshold p-value for significance.
#' @param seed An optional seed.
#' @return A list containing the number of signficant PCs and p value
#' @export
PCPermutationTest = function (dat, B = 100, threshold = 0.05, randomized=F,
                              verbose=TRUE, seed=1, max.pc=100, n.cores=1,
                              center=T, scale=T) {
  ptm = proc.time()
  if(B %% n.cores != 0) stop("Permutations must be an integer multiple of n.cores")
  cat(sprintf("Scaling input matrix [center=%s, scale=%s]\n", center, scale))
  dat = as.matrix(t(scale(t(dat), center=center, scale=scale)))
  set.seed(seed)
  n = min(max.pc, ncol(dat))
  m = nrow(dat)
  print(paste0("Considering only the top ", n, " PCs. Supply max.pc if you wish to change"))
  cat(sprintf("Running initial PCA\n"))
  if(randomized) uu = rsvd::rsvd(dat, k=max.pc)
  else uu = corpcor::fast.svd(dat, tol = 0)
  ndf = n - 1
  dstat = uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  dstat0 = matrix(0, nrow = B, ncol = ndf)
  if(verbose==TRUE) message("Estimating number of significant principal components. Permutation: ")
  if(n.cores==1) {
    for (i in 1:B) {
      if(verbose==TRUE) cat(paste(i, " "))
      dat0 = t(apply(dat, 1, sample, replace = FALSE))
      if(randomized){
        uu0 = rsvd::rsvd(as.matrix(dat0), k=max.pc)
      } else {
        uu0 = corpcor::fast.svd(dat0, tol = 0)
      }
      dstat0[i, ] = uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
  } else {
    cl = parallel::makePSOCKcluster(n.cores, outfile="")
    doParallel::registerDoParallel(cl, n.cores)
    chunksize = B/n.cores
    vals = split(1:B, ceiling(seq_along(1:B)/chunksize))
    dstat0 = foreach::foreach(run.id=1:n.cores, .packages="corpcor", .combine=cbind) %dopar% {
      v = vals[[run.id]]
      #cat(sprintf("Core %s will run perms: %s \n", run.id, paste(v, collapse=",")))
      do.call(rbind, lapply(v, function(i) {
        if(verbose==TRUE) cat(paste(i," "))
        dat0 <- t(apply(dat, 1, sample, replace = FALSE))
        
        if(randomized) uu0 = rsvd::rsvd(as.matrix(dat0), k=max.pc)
        else uu0 = corpcor::fast.svd(dat0, tol = 0)
        uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
      }))
    }
    cat("\nUnregistering parallel backend..")
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    cat(" done\n")
  }
  p = rep(1, n)
  for (i in 1:ndf) {
    p[i] = mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:ndf) {
    p[i] = max(p[(i - 1)], p[i])
  }
  r = sum(p <= threshold)
  y = proc.time() - ptm
  cat(sprintf("\n\n PC permutation test completed. \n %s PCS significant (p<%s, %s bootstraps)\n Runtime: %s s\n ",
              r, threshold, B, signif(y[["elapsed"]], 3)))
  return(list(r = r, p = p))
}

