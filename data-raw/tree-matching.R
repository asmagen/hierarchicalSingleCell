# implement Jiang et al (1994) tree traversal
T_1 <- ReSET::binary_tree_node('a', 
                  left = tree_node('e', left = tree_node('b'), right = tree_node('c')),
                  right = tree_node('d'))

T_2 <- ReSET::binary_tree_node('a', 
                 left = tree_node('b'),
                 right = tree_node('f', left = tree_node('c'), right = tree_node('d')))

# Create the cost matrix
ordered_T1 <- ReSET::postorder(T_1)
ordered_T2 <- ReSET::postorder(T_2)
cost_matrix <- matrix(
  apply(expand.grid(T1 = ordered_T1, T2 = ordered_T2), 1, 
        function(x) abs(as.numeric(charToRaw(x[1])) - as.numeric(charToRaw(x[2])))),
  nrow = length(ordered_T1), byrow = F)
rownames(cost_matrix) <- ordered_T1
colnames(cost_matrix) <- ordered_T2
cost_matrix <- rbind(lambda = rep(0.5, ncol(cost_matrix)), cost_matrix)
cost_matrix <- cbind(lambda = c(0, rep(0.5, nrow(cost_matrix)-1)), cost_matrix)
cost_matrix

initialized <- ReSET::initialize(T_1, T_2, cost_matrix)
