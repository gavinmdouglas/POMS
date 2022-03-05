
filter_rare_table_cols <- function(in_tab, min_nonzero_count, min_nonzero_prop, verbose=TRUE) {
  nonzero_counts <- colSums(in_tab > 0)
  col2remove <- which((nonzero_counts < min_nonzero_count) | (nonzero_counts / nrow(in_tab) < min_nonzero_prop))
  if(length(col2remove) > 0) {
    if(verbose) { message(paste("Filtering", as.character(length(col2remove)), "rare functions from input function table.")) }
    in_tab <- in_tab[, -col2remove]
  }
  return(in_tab)
}


prep_tree_sig_nodes <- function(in_list, taxa_table) {
  
  all_sig_nodes <- c()
  for(func in names(in_list$out_list)) {
    all_sig_nodes <- c(all_sig_nodes,
                       in_list$out_list[[func]]$positive_nodes,
                       in_list$out_list[[func]]$negative_nodes)
  }
  
  all_sig_nodes <- all_sig_nodes[-which(duplicated(all_sig_nodes))]
  
  sig_node_taxa <- list()
  
  for(node in all_sig_nodes) {
    sig_node_taxa[[node]] <- node_taxa(lhs_features = in_list$balances_info$features[[node]]$lhs,
                                               rhs_features = in_list$balances_info$features[[node]]$rhs,
                                               taxa = taxa_table)
  }
  
  tree_sig_subset <- in_list$tree
  nodes2ignore <- which(! tree_sig_subset$node.label %in% all_sig_nodes)
  for(node in all_sig_nodes) {
    tree_sig_subset$node.label[which(tree_sig_subset$node.label == node)] <- sig_node_taxa[[node]]
  }
  tree_sig_subset$node.label[nodes2ignore] <- ""
  
  all_sig_nodes_i <- which(tree_sig_subset$node.label != "") + length(tree_sig_subset$tip.label)
  
  return(list(prepped_tree=tree_sig_subset, nodes2plot=all_sig_nodes_i))
}


prep_tree <- function(phy, tips2keep) {
  # Function to subset tree to tips of MAGs only
  # and run sanity checks.
  # Will midpoint root tree if necessary.
  # Add node labels to tree before returning.
  
  tips2remove <- phy$tip.label[which(! phy$tip.label %in% tips2keep)]
  
  phy <- ape::drop.tip(phy = phy, tip = tips2remove)
  
  if(! ape::is.binary(phy)) {
    stop("Tree is non-binary.")
  }
  
  if(! ape::is.rooted(phy)) {
    phy <- phangorn::midpoint(phy)
  }
  
  phy$node.label <- NULL
  phy <- ape::makeNodeLabel(phy, method="number", prefix='n')
  
  return(phy)
}


subset_abun_table <- function(in_abun, col2keep) {
  
  in_abun <- in_abun[, which(colnames(in_abun) %in% col2keep)]
  
  missing_rows <- which(rowSums(in_abun) == 0)
  missing_samples <- which(colSums(in_abun) == 0)
  
  if(length(missing_rows) > 0) {
    in_abun <- in_abun[-missing_rows, ]
  }
  
  if(length(missing_samples) > 0) {
    in_abun <- in_abun[, -missing_samples]
  }
  
  return(in_abun)
}

