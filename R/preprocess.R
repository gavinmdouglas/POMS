#' @export
filter_rare_table_cols <- function(in_tab, min_nonzero_count, min_nonzero_prop, verbose=TRUE) {
  nonzero_counts <- colSums(in_tab > 0)
  col2remove <- which((nonzero_counts < min_nonzero_count) | (nonzero_counts / ncol(in_tab) < min_nonzero_prop))
  if(length(col2remove) > 0) {
    if(verbose) { message(paste("Filtering", as.character(length(col2remove)), "rare functions from input function table.")) }
    in_tab <- in_tab[, -col2remove]
  }
  return(in_tab)
}

#' @export
prep_tree_sig_nodes <- function(in_list, taxa_table) {
  
  all_sig_nodes <- c()
  for(func in names(in_list$out_list$nodes)) {
    all_sig_nodes <- c(all_sig_nodes,
                       in_list$out_list$nodes[[func]]$positive_nodes,
                       in_list$out_list$nodes[[func]]$negative_nodes)
  }
  
  all_sig_nodes <- all_sig_nodes[-which(duplicated(all_sig_nodes))]
  
  sig_node_taxa <- list()
  
  for(node in all_sig_nodes) {
    sig_node_taxa[[node]] <- node_taxa(lhs_features = in_list$nodes_info$features[[node]]$lhs,
                                               rhs_features = in_list$nodes_info$features[[node]]$rhs,
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
prep_tree_nodes_func <- function(in_list, focal_func, taxa_table) {
  
  all_sig_nodes <- c()
  for(func in names(in_list$out_list$nodes)) {
    all_sig_nodes <- c(all_sig_nodes,
                          in_list$out_list$nodes[[func]]$positive_nodes,
                          in_list$out_list$nodes[[func]]$negative_nodes)
  }
  
  all_sig_nodes <- all_sig_nodes[-which(duplicated(all_sig_nodes))]
  
  sig_node_taxa <- list()
  
  for(node in all_sig_nodes) {
    sig_node_taxa[[node]] <- node_taxa_name(lhs_features = in_list$nodes_info$features[[node]]$lhs,
                                               rhs_features = in_list$nodes_info$features[[node]]$rhs,
                                               taxa = taxa_table)
  }
  
  tree_sig_subset <- in_list$tree
  nodes2ignore <- which(! tree_sig_subset$node.label %in% all_sig_nodes)
  for(node in all_sig_nodes) {
    tree_sig_subset$node.label[which(tree_sig_subset$node.label == node)] <- sig_node_taxa[[node]]
  }
  tree_sig_subset$node.label[nodes2ignore] <- ""
  
  enriched_nodes <- c()
  depleted_nodes <- c()
  
  for(node in all_sig_nodes) {
    
    if(! focal_func %in% rownames(in_list$funcs_per_node[[node]])) {
      next
    }
    
    if(in_list$funcs_per_node[[node]][focal_func, "P"] < 0.05) {
      if(in_list$funcs_per_node[[node]][focal_func, "OR"] > 1) {
        enriched_nodes <- c(enriched_nodes, node)
      } else if(in_list$funcs_per_node[[node]][focal_func, "OR"] < 1) {
        depleted_nodes <- c(depleted_nodes, node)
      }
    }
  }
  
  enriched_nodes_i <- which(in_list$tree$node.label %in% enriched_nodes) + length(in_list$tree$tip.label)
  
  depleted_nodes_i <- which(in_list$tree$node.label %in% depleted_nodes) + length(in_list$tree$tip.label)
  
  all_sig_nodes_minus_enriched_depleted <- all_sig_nodes[which(! all_sig_nodes %in% enriched_nodes)]
  all_sig_nodes_minus_enriched_depleted <- all_sig_nodes_minus_enriched_depleted[which(! all_sig_nodes_minus_enriched_depleted %in% depleted_nodes)]
  
  nonsig_nodes_i <- which(in_list$tree$node.label %in% all_sig_nodes_minus_enriched_depleted) + length(in_list$tree$tip.label)
  
  return(list(prepped_tree=tree_sig_subset, enriched=enriched_nodes_i, depleted=depleted_nodes_i, nonsig=nonsig_nodes_i))
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

#' @export
calc_func_abun <- function(in_abun, in_func, ncores=1) {
  
  out_df <- data.frame(matrix(NA, nrow=ncol(in_func), ncol=ncol(in_abun)))
  colnames(out_df) <- colnames(in_abun)
  rownames(out_df) <- colnames(in_func)
  
  # Check that all rows are found in function table.
  if(length(which(! rownames(in_abun) %in% rownames(in_func))) > 0) {
    stop("Stoppings - some rows in abundance table not found in function table.")
  }
  
  in_func <- in_func[rownames(in_abun), ]
  
  out_sample_func_abun <- mclapply(colnames(in_abun), function(x) { return(colSums(in_abun[, x] * in_func)) }, mc.cores=ncores)
  names(out_sample_func_abun) <- colnames(in_abun)
  
  for(sample in colnames(in_abun)) {
    out_df[, sample] <- out_sample_func_abun[[sample]]
  }
  
  return(out_df)
  
}
