#' @export
filter_rare_table_cols <- function(in_tab, min_nonzero_count, min_nonzero_prop, verbose=TRUE) {
  nonzero_counts <- colSums(in_tab > 0)
  col2remove <- which((nonzero_counts < min_nonzero_count) | (nonzero_counts / nrow(in_tab) < min_nonzero_prop))
  if(length(col2remove) > 0) {
    if(verbose) { message(paste("Filtering", as.character(length(col2remove)), "rare functions from input function table.")) }
    in_tab <- in_tab[, -col2remove]
  }
  return(in_tab)
}

#' @export
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

#' @export
prep_tree_nodes_func <- function(in_list, focal_func, taxa_table = NULL) {
  
  enriched_nodes <- c()
  depleted_nodes <- c()
  nonenriched_sig_nodes <- c()
  nonenriched_nonsig_nodes <- c()
  enriched_nonsig_nodes <- c()
  
  if (!"out_list" %in% names(in_list)) { stop("Stopping - POMS output based on detailed=TRUE is needed.") }
  
  if (!is.null(in_list$out_list[[focal_func]]$positive_nodes)) {
    enriched_nodes <- in_list$out_list[[focal_func]]$positive_nodes
  }
  
  if (!is.null(in_list$out_list[[focal_func]]$negative_nodes)) {
    depleted_nodes <- in_list$out_list[[focal_func]]$negative_nodes
  }
  
  if (!is.null(in_list$out_list[[focal_func]]$nonenriched_sig_nodes)) {
    nonenriched_sig_nodes <- c(nonenriched_sig_nodes, in_list$out_list[[focal_func]]$nonenriched_sig_nodes)
  }
  
  if (!is.null(in_list$out_list[[focal_func]]$missing_sig_nodes)) {
    nonenriched_sig_nodes <- c(nonenriched_sig_nodes, in_list$out_list[[focal_func]]$missing_sig_nodes)
  }
  
  if (!is.null(in_list$out_list[[focal_func]]$missing_nonsig_nodes)) {
    nonenriched_nonsig_nodes <- c(nonenriched_nonsig_nodes, in_list$out_list[[focal_func]]$missing_nonsig_nodes)
  }
  
  if (!is.null(in_list$out_list[[focal_func]]$nonenriched_nonsig_nodes)) {
    nonenriched_nonsig_nodes <- c(nonenriched_nonsig_nodes, in_list$out_list[[focal_func]]$nonenriched_nonsig_nodes)
  }
  
  if (!is.null(in_list$out_list[[focal_func]]$enriched_nonsig_nodes)) {
    enriched_nonsig_nodes <- c(enriched_nonsig_nodes, in_list$out_list[[focal_func]]$enriched_nonsig_nodes)
  }
  
  enriched_nodes_i <- which(in_list$tree$node.label %in% enriched_nodes) + length(in_list$tree$tip.label)
  depleted_nodes_i <- which(in_list$tree$node.label %in% depleted_nodes) + length(in_list$tree$tip.label)
  nonenriched_sig_nodes_i <- which(in_list$tree$node.label %in% nonenriched_sig_nodes) + length(in_list$tree$tip.label)
  nonenriched_nonsig_nodes_i <- which(in_list$tree$node.label %in% nonenriched_nonsig_nodes) + length(in_list$tree$tip.label)
  enriched_nonsig_nodes_i <- which(in_list$tree$node.label %in% enriched_nonsig_nodes) + length(in_list$tree$tip.label)
  sig_nodes_i <- c(enriched_nodes_i, depleted_nodes_i, nonenriched_sig_nodes_i)
  
  in_list$balance_comparisons$mean_direction <- in_list$balance_comparisons$mean_direction[c(in_list$sig_nodes)]
  
  group1_balance_dir <- names(in_list$balance_comparisons$mean_direction)[which(in_list$balance_comparisons$mean_direction == "group1")]
  group1_balance_dir_i <- which(in_list$tree$node.label %in% group1_balance_dir) + length(in_list$tree$tip.label)
  
  group2_balance_dir <- names(in_list$balance_comparisons$mean_direction)[which(in_list$balance_comparisons$mean_direction == "group2")]
  group2_balance_dir_i <- which(in_list$tree$node.label %in% group2_balance_dir) + length(in_list$tree$tip.label)
  
  lhs_balance_function_enrich_dir_i <- c(group1_balance_dir_i[which(group1_balance_dir_i %in% enriched_nodes_i)],
                                         group2_balance_dir_i[which(group2_balance_dir_i %in% depleted_nodes_i)])
  
  rhs_balance_function_enrich_dir_i <- c(group2_balance_dir_i[which(group2_balance_dir_i %in% enriched_nodes_i)],
                                         group1_balance_dir_i[which(group1_balance_dir_i %in% depleted_nodes_i)])
  
  tested_nodes_i <- which(in_list$tree$node.label %in% names(in_list$balances_info$balances)) + length(in_list$tree$tip.label)
  
  # Specify taxa names at significant nodes if taxa table provided.
  if (!is.null(taxa_table)) {
    tested_nodes <- names(in_list$balances_info$balances)
    sig_node_taxa <- list()
    for (node in tested_nodes) {
      sig_node_taxa[[node]] <- node_taxa(lhs_features = in_list$balances_info$features[[node]]$lhs,
                                         rhs_features = in_list$balances_info$features[[node]]$rhs,
                                         taxa = taxa_table)
    }
    tree_sig_subset <- in_list$tree
    nodes2ignore <- which(!tree_sig_subset$node.label %in% tested_nodes)
    for (node in tested_nodes) {
      tree_sig_subset$node.label[which(tree_sig_subset$node.label == node)] <- sig_node_taxa[[node]]
    }
    tree_sig_subset$node.label[nodes2ignore] <- ""
  } else {
    tree_sig_subset <- in_list$tree
  }
  
  return(list(prepped_tree = tree_sig_subset,
              group1_func_enriched_nodes = enriched_nodes_i,
              group2_func_enriched_nodes = depleted_nodes_i,
              nonenriched_sig = nonenriched_sig_nodes_i,
              sig_balances_nodes = sig_nodes_i,
              nonenriched_nonsig = nonenriched_nonsig_nodes_i,
              enriched_nonsig = enriched_nonsig_nodes_i,
              tested_nodes = tested_nodes_i,
              group1_balance_dir_i = group1_balance_dir_i,
              group2_balance_dir_i = group2_balance_dir_i,
              lhs_balance_function_enrich_dir_i = lhs_balance_function_enrich_dir_i,
              rhs_balance_function_enrich_dir_i = rhs_balance_function_enrich_dir_i))
}



#' @export
prep_tree <- function(phy, tips2keep) {
  # Function to subset tree to tips of MAGs only
  # and run sanity checks.
  # Will midpoint root tree if necessary.
  # Add node labels to tree before returning.
  
  tips2remove <- phy$tip.label[which(! phy$tip.label %in% tips2keep)]
  
  phy <- ape::drop.tip(phy = phy, tip = tips2remove)
  
  if(! ape::is.binary.tree(phy)) {
    stop("Tree is non-binary.")
  }
  
  if(! ape::is.rooted(phy)) {
    phy <- phangorn::midpoint(phy)
  }
  
  phy$node.label <- NULL
  phy <- ape::makeNodeLabel(phy, method="number", prefix='n')
  
  return(phy)
}

#' @export
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
