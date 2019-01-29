T_1 <- ReSET::tree_node('a', children = list(ReSET::tree_node('b'),
                                             ReSET::tree_node('c'), 
                                             ReSET::tree_node('d')))

T_2 <- ReSET::tree_node('a', children = list(ReSET::tree_node("b"), 
                                             ReSET::tree_node("c")))
par(mfrow = c(1, 2))
plot(T_1)
plot(T_2)

ordered_T1 <- ReSET::postorder_labels2(T_1)
ordered_T2 <- ReSET::postorder_labels2(T_2)
cost_matrix <- matrix(apply(expand.grid(T1 = ordered_T1, T2 = ordered_T2), 1, 
                            function(x) abs(as.numeric(charToRaw(x[1])) - as.numeric(charToRaw(x[2])))),
                      nrow = length(ordered_T1), byrow = F)
rownames(cost_matrix) <- ordered_T1
colnames(cost_matrix) <- ordered_T2
cost_matrix <- rbind(lambda = rep(0.3, ncol(cost_matrix)), cost_matrix)
cost_matrix <- cbind(lambda = c(0, rep(0.3, nrow(cost_matrix)-1)), cost_matrix)
cost_matrix

align_obj <- ReSET::create_align_object2(T_1, T_2, cost_matrix)
align_obj <- ReSET::fill_matrix2(align_obj)

align_obj <- ReSET::align(T_1, T_2, cost_matrix)
plot(align_obj$tree)
