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
#' @param FC_thresh a positive real number that defines the cutoff of fold change the DE analysis
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

#' a wrapper function of DE_childrenNodes that finds the overlapping DE genes
#'
#' @param X_exp1, X_exp2 n genes by m cell count matrix of gene expression for two experiments
#' @param cn_1_idx_exp1, cn_2_idx_exp2 a m by 1 logical vector that indicates where the cells for children node 1 are for experiement 1 and 2
#' @param cn_2_idx_exp1, cn_2_idx_exp2 a m by 1 logical vector that indicates where the cells for children node 2 are for experiment 1 and 2
#' @param q_thresh_1, q_thresh_2 a positive real number that gives cutoff of q value of DE analysis for experiment 1 and 2
#' @param FC_thresh_1, FC_thresh_2 a positive real number that defines the cutoff of fold change the DE analysis for experiment 1 and 2
#'
#' @return the overlapping differentially expressed genes between experiment 1 and 2

DE_twoExps = function(X_exp1, cn_1_idx_exp1, cn_2_idx_exp1, q_thresh_exp1, FC_thresh_exp1,
                      X_exp2, cn_1_idx_exp2, cn_2_idx_exp2, q_thresh_exp2, FC_thresh_exp2) {
       DE_hits_exp1 = DE_childrenNodes(X_exp1, cn_1_idx_exp1, cn_2_idx_exp1, q_thresh_exp1, FC_thresh_exp1)
       DE_hits_exp2 = DE_childrenNodes(X_exp2, cn_1_idx_exp2, cn_2_idx_exp2, q_thresh_exp2, FC_thresh_exp2)
       colnames(DE_hits_exp1) = c('log2FoldChange_exp1', 'padj_exp1')
       colnames(DE_hits_exp2) = c('log2FoldChange_exp2', 'padj_exp2')
       overlap_hits = intersect(rownames(DE_hits_exp1), rownames(DE_hits_exp2))
       DE_hits_overlap = cbind(DE_hits_exp1[overlap_hits,], DE_hits_exp2[overlap_hits,])
       return(DE_hits_overlap)
}