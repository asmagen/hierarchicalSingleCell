suppressMessages({
    library(tidyverse)
    library(DESeq2)
})

GeneSignature = function(data, label, alpha=0.05, out.path=NULL) {
    dds <- DESeqDataSetFromMatrix(countData = data, # genes x samples
                              colData = data.frame(is_community = label),
                              design= ~ is_community)
    dds <- DESeq(dds)
    resultsNames(dds)
    res <- results(dds, alpha=alpha)
    resOrdered <- res[order(res$pvalue),]
    sig.idx = res$padj < alpha
    print(paste0(sum(sig.idx), 'genes found differentially expressed')
    if(!is.null(out.path)) write.csv(res[sig.idx,], file=out.path)
}

data_3k_path = '/home/ubuntu/robustSingleCell/data/t_3k_genesymbol.csv'
#data_4k_path = '/home/ubuntu/robustSingleCell/data/t_4k_genesymbol.csv'
data = read.table(data_3k_path)  
gene_names = data[,1]
data = data[,-1]
rownames(data) = gene_names
         
# make up some random labels
n_communities = 5
n_cells = dim(data)[2]
communities = sample(seq(n_communities), n_cells, replace=TRUE)
          
# run the 
base.path = '/home/ubuntu/robustSingleCell/data/DESeq'
for (i in unique(communities)) {
    suppressMessages(GeneSignature(t(test_data), label=(communities == i) %>% as.numeric, out.path=paste0(base.path, 'community_', i, '.csv')))
}