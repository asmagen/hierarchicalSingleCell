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
