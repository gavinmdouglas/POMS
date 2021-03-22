#' @export
calc_node_and_func_stats <- function(POMS_out, group1_samples, group2_samples) {
  
  pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(POMS_out$balances_info$balances,
                                                            group1_samples,
                                                            group2_samples,
                                                            skip_wilcoxon=FALSE)
  
  node_effect_size <- sapply(POMS_out$balances_info$balances,
                             function(x) {
                               mean_diff <- mean(x[group1_samples]) - mean(x[group2_samples])
                               pooled_sd <- sqrt((sd(x[group1_samples]) ** 2 + sd(x[group2_samples]) ** 2) / 2)
                               return(mean_diff / pooled_sd)
                             })
  
  func_OR_table <- data.frame(matrix(NA,
                                     nrow=length(POMS_out$balances_info$balances),
                                     ncol=nrow(POMS_out$df)))
  colnames(func_OR_table) <- rownames(POMS_out$df)
  rownames(func_OR_table) <- names(POMS_out$balances_info$balances)
  
  for(node in names(POMS_out$funcs_per_node)) {
    for(func in rownames(POMS_out$funcs_per_node[[node]])) {
      if(pairwise_node_out$mean_direction[node] == "group2") {
        func_OR_table[node, func] <- 1 / POMS_out$funcs_per_node[[node]][func, "OR"]
      } else {
        func_OR_table[node, func] <- POMS_out$funcs_per_node[[node]][func, "OR"]
      }
      
    }
  }
  
  func_P_table <- data.frame(matrix(NA,
                                    nrow=length(POMS_out$balances_info$balances),
                                    ncol=nrow(POMS_out$df)))
  colnames(func_P_table) <- rownames(POMS_out$df)
  rownames(func_P_table) <- names(POMS_out$balances_info$balances)
  
  for(node in names(POMS_out$funcs_per_node)) {
    for(func in rownames(POMS_out$funcs_per_node[[node]])) {
      func_P_table[node, func] <- POMS_out$funcs_per_node[[node]][func, "P"]
    }
  }
  
  return(list(node_wilcox_p=pairwise_node_out$wilcox_raw_p,
              node_effect_size=node_effect_size,
              func_OR_table=func_OR_table,
              func_P_table=func_P_table))
}

#' @export
filter_out_nonsig_funcs <- function(in_stats, min_sig_nodes=5) {
  
  funcs2remove <- which(colSums(in_stats$func_P_table < 0.05, na.rm = TRUE) < min_sig_nodes)
  
  if(length(funcs2remove == ncol(in_stats))) {
    stop("All functions removed after filtering functions not significantly different at specified number of nodes.")
  } else if(length(funcs2remove) > 0) {
    in_stats$func_OR_table <- in_stats$func_OR_table[ , -funcs2remove, drop=FALSE]
    in_stats$func_P_table <- in_stats$func_P_table[ , -funcs2remove, drop=FALSE]
  }
  
  return(in_stats)
}


#' @export
POMS_multinomial_test <- function(in_stats, node_p_cutoff=0.05, gene_p_cutoff=0.05) {
  
  num_tested_nodes <- length(in_stats$node_wilcox_p)
  
  prop_sig_node_balances <- length(which(in_stats$node_wilcox_p < node_p_cutoff)) / num_tested_nodes
  
  sig_nodes <- names(in_stats$node_wilcox_p)[which(in_stats$node_wilcox_p < node_p_cutoff)]
  sig_nodes_pos <- names(which(in_stats$node_effect_size[sig_nodes] > 0))
  sig_nodes_neg <- names(which(in_stats$node_effect_size[sig_nodes] < 0))
  
  multinom_results <- data.frame(matrix(NA, nrow=ncol(in_stats$func_P_table), ncol=0))
  rownames(multinom_results) <- colnames(in_stats$func_P_table)
  
  for(gene in colnames(in_stats$func_P_table)) {
    
    sig_nodes_for_gene <- rownames(in_stats$func_P_table[which(in_stats$func_P_table[, gene] < gene_p_cutoff), , drop=FALSE])
    
    sig_nodes_func_OR <- in_stats$func_OR_table[sig_nodes_for_gene, gene, drop=TRUE]
    
    pos_sig_nodes_for_gene <- sig_nodes_for_gene[which(sig_nodes_func_OR > 1)]
    
    neg_sig_nodes_for_gene <- sig_nodes_for_gene[which(sig_nodes_func_OR < 1)]
    
    sig_same_dir_nodes <- c(pos_sig_nodes_for_gene[which(pos_sig_nodes_for_gene %in% sig_nodes_pos)],
                            neg_sig_nodes_for_gene[which(neg_sig_nodes_for_gene %in% sig_nodes_neg)])
    
    sig_diff_dir_nodes <- c(pos_sig_nodes_for_gene[which(pos_sig_nodes_for_gene %in% sig_nodes_neg)], neg_sig_nodes_for_gene[which(neg_sig_nodes_for_gene %in% sig_nodes_pos)])
    
    num_sig_same_dir_nodes <- length(sig_same_dir_nodes)
    
    num_sig_diff_dir_nodes <- length(sig_diff_dir_nodes)
    
    # Note that the expected prop are the same for the same / diff because they are assumed to be 50/50.
    expected_prop <- c(prop_sig_node_balances * 0.5, prop_sig_node_balances * 0.5, 1 - prop_sig_node_balances)
    observed_counts <- c(num_sig_same_dir_nodes, num_sig_diff_dir_nodes,
                         length(sig_nodes_for_gene) - num_sig_same_dir_nodes - num_sig_diff_dir_nodes)
    
    multinom_results[gene, "num_sig_nodes_balances"] <- length(sig_nodes)
    multinom_results[gene, "num_sig_nodes_genes"] <- length(sig_nodes_for_gene)
    multinom_results[gene, "num_tested_nodes"] <- num_tested_nodes
    multinom_results[gene, "exp_same.intersect_prop"] <- prop_sig_node_balances * 0.5
    multinom_results[gene, "exp_diff.intersect_prop"] <- prop_sig_node_balances * 0.5
    multinom_results[gene, "exp_nonintersect_prop"] <- 1 - prop_sig_node_balances
    multinom_results[gene, "obs_same.intersect_count"] <- num_sig_same_dir_nodes
    multinom_results[gene, "obs_diff.intersect_count"] <- num_sig_diff_dir_nodes
    multinom_results[gene, "obs_nonintersect_count"] <- length(sig_nodes_for_gene) - num_sig_same_dir_nodes - num_sig_diff_dir_nodes
    
    multinom_results[gene, "multinom_p"] <- xmulti(obs=observed_counts,
                                                   expr=expected_prop, detail=0)$pProb
    
  }
  
  multinom_results$multinom_BH <- p.adjust(multinom_results$multinom_p, "BH")
  
  return(multinom_results)
}