#' Given a ancestor matrix of the hierarchy and corresponding description
#' of the leaves, return a name for each internal nodes
#' @param id_matrix a m by k matrix that indicates the membership of
#' nodes at each level ($m$ is the number of cells and $k$ is the number of 
#' splits)
#' @param cell_desc a vector of size m that indicates the annotation at the
#' bottom of the hierarchy
#' @return a ancestor matrix replacing id with corresponding name according to
#' last number
word_ancestor_matrix = function(id_matrix, cell_desc){
    K = seq(ncol(id_matrix))
    desc.df = c()
    for (k in K) {
        id = id_matrix[,k]
        unique_labels = data.frame(id = id, desc=cell_desc) %>% distinct()
        label.df = concatLabels(unique_labels)
        label.map = hashmap(label.df$id, label.df$desc)
        desc = sapply(id, function(x) label.map[[as.numeric(x)]])
        desc.df = cbind(desc.df, desc)
    }                
    colnames(desc.df) = K
    return(desc.df)
}

#' A helper function that takes node id and return the name of all leaves that
#' are its children
#' @param unique_labels a c by 2 data frame with column names id (cluster id) and
#' desc (the name of any leave that is its child)
#' @return a data frame that contains node id and concatenated names of all its children
#' leaves
concatLabels = function(unique_labels) {
    label.df = unique_labels %>%
      group_by(id) %>%
      mutate(n = seq(n())) %>%
      spread(n, desc) %>%
      unite(desc, colnames(.)[-1]) %>%
      mutate(desc = sapply(desc, function(x) gsub('_NA', '', x)))
    return(label.df)
}


