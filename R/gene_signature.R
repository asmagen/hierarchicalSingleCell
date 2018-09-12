suppressMessages({
    library(tidyverse)
    library(DESeq2)
    library(limma)
})

#' Return the differentially expressed genes between two children node of a tree
#'
#' @param n genes by m cell count matrix of gene expression
#' @param child_node_1_idx a m by 1 logical vector that indicates where the cells for children node 1 are
#' @param child_node_2_idx a m by 1 logical vector that indicates where the cells for children node 2 are
#' @param q_thresh a positive real number that gives cutoff of q value of DE analysis
#' @param FC_thresh a positive real number that defines the cutoff of fold change the DE analysis
#' 
#' @return a k by 2 matrix that contains k differentially expressed genes between
#' the two child nodes of corresponding q value and FC value with rownames of the
#' genes
#' Usage:
#' X = matrix(sample(1:100, 5000, replace=T), nrow=100) # 100 genes x 50 cells
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
                                     colData = data.frame(children=factor(children_label)),
                                     design = ~children)
    suppressMessages({DE_res = DESeq(DE_data, parallel=TRUE)})
    DE_hits = results(DE_res, alpha=q_thresh)
    DE_hits = na.omit(DE_hits[order(DE_hits$log2FoldChange, decreasing=T),])
    hits_idx = DE_hits$log2FoldChange > log(FC_thresh, base=2)
    DE_hits_final = DE_hits[hits_idx,c(2,6)]
    rownames(DE_hits_final) = rownames(DE_hits)[hits_idx]
    return(data.matrix(DE_hits_final))
}

DE_childrenNodes.limma = function(X, child_node_1_idx, child_node_2_idx,
                                  q_thresh, FC_thresh) {
    X_child_node_1 = X[,child_node_1_idx]
    X_child_node_2 = X[,child_node_2_idx]
    X_children_node = cbind(X_child_node_1, X_child_node_2)
    children_label = c(rep(1, sum(child_node_1_idx)), rep(2, sum(child_node_2_idx)))
    design = model.matrix(~children_label)
    my.lm = limma::lmFit(X, design = design) # X is a normalized data matrix
    my.lm = eBayes(my.lm) 
    FC_genes = topTable(my.lm, number=Inf, coef = seq(2,ncol(design)))
    FC_genes = FC_genes[FC_genes$adj.P.Val < q_thresh & FC_genes$logFC > log(FC_thresh, base=2), c('logFC', 'adj.P.Val')]
    FC_genes = FC_genes[order(FC_genes$logFC, decreasing =T),]
    return(FC_genes)
}

DE_childrenNodes.wilcox = function(X, child_node_1_idx, child_node_2_idx,
                                  q_thresh, FC_thresh) {
    X_child_node_1 = X[,child_node_1_idx]
    X_child_node_2 = X[,child_node_2_idx]
    X_children_node = cbind(X_child_node_1, X_child_node_2)
    children_label = c(rep(1, sum(child_node_1_idx)), rep(2, sum(child_node_2_idx)))
    test_matrix = apply(X, 1, function(x) {
        test.res = wilcox.test(x ~ children_label)
        p.value = test.res$p.value
        logFC = log(mean(x[child_node_1_idx])/mean(x[child_node_2_idx]), base=2)
        c(logFC, p.value)
    })
    test_matrix = t(test_matrix)
    #print(test_matrix)
    q.value = qvalue(test_matrix[,2])$qvalues
    test_matrix = data.frame(cbind(test_matrix, q.value))
    colnames(test_matrix) = c('logFC', 'p.value', 'q.value')
    rownames(test_matrix) = rownames(X)
    test_matrix[q.value < q_thresh & test_matrix[,1] > log(FC_thresh, base=2),]
}

#' a wrapper function of DE_childrenNodes that finds the overlapping DE genes
#'
#' @param X_exp1, X_exp2 n genes by m cell count matrix of gene expression for two experiments
#' @param cn_1_idx_exp1, cn_2_idx_exp2 a m by 1 logical vector that indicates 
#' where the cells for children node 1 are for experiement 1 and 2
#' @param cn_2_idx_exp1, cn_2_idx_exp2 a m by 1 logical vector that indicates 
#' where the cells for children node 2 are for experiment 1 and 2
#' @param q_thresh_1, q_thresh_2 a positive real number that gives cutoff of q 
#' value of DE analysis for experiment 1 and 2
#' @param FC_thresh_1, FC_thresh_2 a positive real number that defines the cutoff 
#' of fold change the DE analysis for experiment 1 and 2
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

#' Wrapper given expression dataset and a logical vector of membership
#' 
#' @param membership a binary vector of cluster assignment
#' @param expression the n genes by m cells that corresponding to the 
#' cells in the cluster membership 
#' @param q_thresh FDR threshold
#' @param FC_thresh FC threshold
#' @return the gene symbol of differentially expressed genes
OverExpressedGenes <- function(membership, expression, q_thresh, FC_thresh) {
    DE_res = DE_childrenNodes.limma(expression, membership, !membership, q_thresh=q_thresh, FC_thresh=FC_thresh)
    return(rownames(DE_res))
}

#' @param tree a m cells by k hierarchy data matrix that indicates the hierarchy 
#' structure
#' @param expression a n genes by m cells data matrix that records raw gene coutns
#' @param q_thresh FDR threshold
#' @param FC_thresh FC threshold
#' @return the set of differentially expressed genes for each intermediate
#' node in the hierarchy
get.node.genes <- function(tree, expression, q_thresh, FC_thresh, k=1) {
    if (ncol(tree) < k) {
        return(NULL)
    } else {
        membership_raw = tree[,k]
        membership = membership_raw == unique(membership_raw)[1]
        over1 = OverExpressedGenes(membership, expression, q_thresh, FC_thresh)        
        subtree1_id = paste(data.matrix(tree[membership, seq(k)])[1,], collapse='_')
        print(paste(subtree1_id, length(over1)))
        if(length(over1) > 0) {
            combined = data.frame(id=subtree1_id, gene=over1)
        } else {
            combined = NULL
        }
        over2 = OverExpressedGenes(!membership, expression, q_thresh, FC_thresh)
        subtree2_id = paste(data.matrix(tree[!membership, seq(k)])[1,], collapse='_')
        print(paste(subtree2_id, length(over2)))
        if(length(over2) > 0) {
            combined = rbind(combined, data.frame(id=subtree2_id, gene=over2))
        } 
        
        combined.left = get.node.genes(data.matrix(tree[membership,]), expression[,membership], q_thresh, FC_thresh, k=k+1)
        if(!is.null(combined.left)) combined = rbind(combined, combined.left)
        combined.right = get.node.genes(data.matrix(tree[!membership,]), expression[,!membership], q_thresh, FC_thresh, k=k+1)
        if(!is.null(combined.right)) combined = rbind(combined, combined.right)

        return(combined)
    }
}

