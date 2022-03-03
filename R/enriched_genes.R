feature_sets_func_fisher <- function(in_func, feature_set1, feature_set2, add_pseudocount=FALSE, multiple_test_corr="none") {
  
  in_func_set1 <- in_func[feature_set1, , drop = FALSE]
  in_func_set2 <- in_func[feature_set2, , drop = FALSE]
  
  func_ids_all <- colnames(in_func)
  
  set1_func_missing <- func_ids_all[which(colSums(in_func_set1[, func_ids_all, drop = FALSE]) == 0)]
  set2_func_missing <- func_ids_all[which(colSums(in_func_set2[, func_ids_all, drop = FALSE]) == 0)]
  both_func_missing <- set1_func_missing[which(set1_func_missing %in% set2_func_missing)]
  
  if(length(both_func_missing) > 0) {
    in_func_set1 <- in_func_set1[, -which(colnames(in_func_set1) %in% both_func_missing), drop = FALSE]
    in_func_set2 <- in_func_set2[, -which(colnames(in_func_set2) %in% both_func_missing), drop = FALSE]
  }
  
  func_ids <- colnames(in_func_set1)
  
  fisher_out <- data.frame(matrix(NA, nrow=length(func_ids), ncol=6))
  colnames(fisher_out) <- c("set1_pos", "set1_neg", "set2_pos", "set2_neg", "OR", "P")
  
  rownames(fisher_out) <- func_ids
  
  for (func in func_ids) {
    
    func_count_matrix <- matrix(c(length(which(in_func_set1[, func] > 0)), length(which(in_func_set1[, func] == 0)),
                                  length(which(in_func_set2[, func] > 0)), length(which(in_func_set2[, func] == 0))),
                                nrow=2, ncol=2)
    
    if (add_pseudocount) {
      func_count_matrix <- func_count_matrix + 1
    }
    
    func_count_fisher <- fisher.test(round(func_count_matrix))
    
    fisher_out[func, ] <- c(func_count_matrix[1, 1], func_count_matrix[2,1], func_count_matrix[1, 2], func_count_matrix[2,2],
                            func_count_fisher$estimate, func_count_fisher$p.value)
  }
  
  fisher_out$P_corr <- p.adjust(fisher_out$P, multiple_test_corr)
  
  return(fisher_out)
}

node_func_fisher <- function(node, in_tree, in_func, higher_group, pseudocount=NULL, multiple_test_corr="none") {
  node_features <- lhs_rhs_asvs(in_tree, node, get_node_index=TRUE)
  if(higher_group == "group1") {
    node_fisher_tests <- feature_sets_func_fisher(in_func = in_func,
                                                  feature_set1 = node_features$lhs,
                                                  feature_set2 = node_features$rhs,
                                                  pseudocount=pseudocount,
                                                  multiple_test_corr=multiple_test_corr)
  } else if(higher_group == "group2") {
    node_fisher_tests <- feature_sets_func_fisher(in_func = in_func,
                                                  feature_set1 = node_features$rhs,
                                                  feature_set2 = node_features$lhs,
                                                  pseudocount=pseudocount,
                                                  multiple_test_corr=multiple_test_corr)
  }
  return(node_fisher_tests)
}

calc_OR <- function(exposed_pos, exposed_neg, control_pos, control_neg) {
  return((exposed_pos / exposed_neg) / (control_pos / control_neg))  
}



