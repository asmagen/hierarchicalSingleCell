# T3K and T4K data
load('inst/extdata/T3k.RData')
membership <- tree[,9]
colnames(normalized) <- NULL
T3K <- list(normalized_data = normalized, membership = membership)
devtools::use_data(T3K)

load('inst/extdata/T4k.RData')
membership <- tree[,11]
membership[membership == 'CD4_Th1.2'] <- 'CD4.Th1.2'
colnames(normalized) <- NULL
T4K <- list(normalized_data = normalized, membership = membership)
devtools::use_data(T4K)
