source('../../R/gene_signature.R')
print('reading in hierarchy...')
t3.hierarchy.data = read.table('../../data/treeHierarchy/t3.scaledexpression.cutree.csv',
                    row.names=1)
CD4.idx = t3.hierarchy.data[,'CD4'] > 0
CD8alpha.idx = t3.hierarchy.data[,'CD8A'] > 0
subset_cells.idx = (CD4.idx | CD8alpha.idx)
t3.hierarchy.data = t3.hierarchy.data[subset_cells.idx,]

hierarchy = t3.hierarchy.data[,c('X2', 'X4', 'X8', 'X16')]
print('getting the normalized expression matrix...')
t3.data = t3.hierarchy.data[,!(colnames(t3.hierarchy.data) %in% c('X2', 'X4', 'X8', 'X16'))]
t3.data.arranged = t(t3.data)

print('running the differentially expressed data across hierarchy..')
res = get.node.genes(hierarchy, t3.data.arranged, 0.1, 1.5)
write.table(combined, file='../../data/treeHierarchy/t3.diff.csv', row.names=F)
