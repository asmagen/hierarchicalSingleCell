suppressMessages({
    library(tidyverse)
    library(DESeq2)
})

#' Return the differentially expressed genes between two children node of a tree
#'
#' @param X n genes by m cell count matrix of gene expression
#' @param child_node_1_idx a m by 1 logical vector that indicates where the cells for children node 1 are
#' @param child_node_2_idx a m by 1 logical vector that indicates where the cells for children node 2 are
#' @param q_thresh a positive real number that gives cutoff of q value of DE analysis
#' @FC_thresh a positive real number that defines the cutoff of fold change the DE analysis
#' 
#' @return a k by 2 matrix that contains k differentially expressed genes between the two child nodes of corresponding q value and FC value with rownames of the genes
#' Usage:
#' X = matrix(sample(1:100, n=5000), nrow=100) # 100 genes x 50 cells
#' child_node_1_idx = sample(c(T, F), 50, replace=T)
#' child_node_2_idx = !child_node_1_idx
#' q_thresh = 0.05
#' FC_thresh = 1
#' DE_childrenNodes(X, child_node_1_idx, child_node_2_idx, q_thresh, FC_thresh)

DE_childrenNodes = function(X, child_node_1_idx, child_node_2_idx, q_thresh, FC_thresh) {
    X_child_node_1 = X[,child_node_1_idx]
    X_child_node_2 = X[,child_node_2_idx]
    X_children_node = cbind(X_child_node_1, X_child_node_2)
    children_label = c(rep(1, sum(child_node_1_idx)), rep(2, sum(child_node_2_idx)))
    DE_data = DESeqDataSetFromMatrix(countData = X_children_node,
                                     colData = data.frame(children=children_label),
                                     design = ~children)
    DE_res = DESeq(DE_data, parallel=TRUE)
    DE_hits = results(DE_res, alpha=q_thresh)
    DE_hits = na.omit(DE_hits[order(DE_hits$log2FoldChange, decreasing=T),])
    DE_hits_final = DE_hits[DE_hits$log2FoldChange > log(FC_thresh, base=2),]
    return(data.matrix(DE_hits_final[,c(2, 6)]))
}