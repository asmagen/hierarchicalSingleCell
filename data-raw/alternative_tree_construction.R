library(clustree)
clusters.aggregate = cbind(1,clusters. cluster.membership.vector)
colnames(clusters.aggregate) = paste('K',seq(ncol(clusters.aggregate)))
print(clustree(clusters.aggregate, prefix = "K"))
