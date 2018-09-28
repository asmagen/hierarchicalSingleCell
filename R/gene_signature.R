#' Differential expression that iteratively compares one cluster against all the others
#' @param bottom_id the cluster id of each cell
#' @param normalized the normalized expression matrix with column ordered as the bottom_id
#' @param q_thresh q value threshold for DE analysis
#' @param FC_thresh fold change threshold for DE analysis
#' @return the list of genes that are overexpressed for each cluster against all others 
SimpleDE = function(bottom_id, normalized, q_thresh = 0.1, FC_thresh = 1.2) {
    FC_genes_all = c()
    bottom_unique_id = unique(bottom_id)
    for (unique_id in bottom_unique_id) {
        is_in_id = bottom_id == unique_id
        design = model.matrix(~is_in_id)
        my.lm = limma::lmFit(normalized, design=design)
        my.lm = eBayes(my.lm)
        FC_genes = topTable(my.lm, number=Inf, coef = seq(2,ncol(design)))
        FC_genes = FC_genes[FC_genes$adj.P.Val < q_thresh, c('logFC', 'adj.P.Val')]
        FC_genes = FC_genes[FC_genes$logFC > log(FC_thresh, base=2),]
        FC_genes = FC_genes[order(FC_genes$logFC, decreasing=T),]
        FC_genes$id = unique_id
        FC_genes_all = rbind(FC_genes_all, FC_genes)
    }
    return(FC_genes_all)
}

#' Use limma to find differentially expressed genes
#' @param expression m by n expression matrix with $m$ genes and $n$ cells
#' @param unique_cluster_membership a logical vector of size $n$ that indicates membership
#' @param q_thresh q value threshold for DE analysis
#' @param FC_thresh fold change threshold for DE analysis
sensDE_limma = function(expression, unique_cluster_membership, q_thresh = 0.1, FC_thresh = 1.2) {
    design = model.matrix(~unique_cluster_membership)
    my.lm = limma::lmFit(expression, design=design)
    my.lm = eBayes(my.lm)

    FC_genes = topTable(my.lm, number=Inf, coef = seq(2,ncol(design)))
    FC_genes = FC_genes[FC_genes$adj.P.Val < q_thresh & FC_genes$logFC > log(FC_thresh, base=2), c('logFC', 'adj.P.Val')]
    FC_genes = FC_genes[order(FC_genes$logFC, decreasing =T),]
    return(FC_genes)
}


#' sensitive DE at level k of the hierarchy
#' @param cluster_name the name of the clusters that each cell belongs to
#' @param expression m by n expression matrix for $m$ genes and $n$ cells
#' @param q_thresh q value threshold for DE analysis
#' @param FC_thresh fold change threshold for DE analysis
#' @return differentially expressed genes for each cluster
DE_level_k = function(cluster_name, expression, q_thresh, FC_thresh) {
    unique_clusters = unique(cluster_name)
    for (unique_cluster in unique_clusters) {
        unique_cluster_membership = (cluster_name == unique_cluster)
        FC_genes = sensDE_limma(expression, unique_cluster_membership, q_thresh, FC_thresh)
        FC_genes$node = unique_cluster
        cluster_DE = rbind(cluster_DE, FC_genes)
    }
    return(cluster_DE)
}

#' sensitive DE at all levels of the hierarchy
#' @param expression.all the $m$ genes by $n$ cells expression matrix
#' @param cell_ancestor_matrix hierarchy of the cells
#' @param q_thresh q value threshold for DE analysis
#' @param FC_thresh fold change threshold for DE analysis
#' @return differentially expressed genes for each node in the hierarchy
sens_DE = function(expression.all, cell_ancestor_matrix, q_thresh, FC_thresh) {
    cluster_DE_all = c()
    for (k in K) {
    if (k == 1) {
        expression = expression.all
        cluster_name = cell_ancestor_matrix[,k]
    } else {
        summary_matrix = table(cell_ancestor_matrix[,k-1], cell_ancestor_matrix[,k])
        parent_split_idx = apply(summary_matrix, 1, function(x) sum(x > 0) == 2)
        parent_split = names(parent_split_idx)[parent_split_idx]
        subset_idx = cell_ancestor_matrix[, k-1] == parent_split
        expression = expression.all[, subset_idx]
        cluster_name = cell_ancestor_matrix[subset_idx, k]
    }
    cluster_DE_k = DE_level_k(cluster_name, expression, q_thresh, FC_thresh)
    cluster_DE_all = rbind(cluster_DE_all, cluster_DE_k)
    }
    return(cluster_DE_all)
}
