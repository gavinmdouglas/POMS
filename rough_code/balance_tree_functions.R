library(ape)
library(parallel)
library(phangorn)
library(stringr)

two_group_balance_tree_pipeline <- function(abun,
                                            func,
                                            phylogeny,
                                            taxa,
                                            group1_samples,
                                            group2_samples,
                                            ncores=1,
                                            pseudocount=1,
                                            significant_nodes=NULL,
                                            tested_balances=NULL,
                                            min_num_tips=10,
                                            min_func_instances=10,
                                            min_func_prop=0.001,
                                            name_threshold=0.9,
                                            balance_p_cutoff = 0.05,
                                            balance_correction = "BY",
                                            function_p_cutoff = 0.05,
                                            function_correction = "none",
                                            func_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_descrip_22Aug2019.tsv",
                                            verbose=FALSE) {
  
  if(verbose) { message("Checking input arguments.") }
  input_param <- check_two_group_balance_args(abun=abun, func=func, phylogeny=phylogeny, taxa=taxa,
                                              group1_samples=group1_samples, group2_samples=group2_samples,
                                              ncores=ncores, pseudocount=pseudocount, significant_nodes=significant_nodes,
                                              tested_balances=tested_balances, min_num_tips=min_num_tips, min_func_instances=min_func_instances,
                                              min_func_prop=min_func_prop, name_threshold=name_threshold,
                                              balance_p_cutoff=balance_p_cutoff, balance_correction=balance_correction,
                                              function_p_cutoff=function_p_cutoff, function_correction=function_correction,
                                              func_descrip_infile=func_descrip_infile, verbose=verbose)
  
  if(verbose) { message("Prepping input phylogeny.") }
  phylogeny <- prep_tree(phy=phylogeny, tips2keep=rownames(abun))
  
  if((min_func_instances > 0) || (min_func_prop > 0)) {
    func <- filter_rare_table_cols(func,  min_func_instances, min_func_prop, verbose)
  }
  
  if(is.null(significant_nodes)) {
    
    if(verbose) { message("Calculating balances.") }
    
    calculated_balances <- compute_tree_node_balances(abun=abun, phylogeny=phylogeny, ncores=ncores, min_num_tips = min_num_tips)
    
    if(verbose) { message("Identified ", length(calculated_balances$balances), " of ", length(phylogeny$node.label), " nodes as non-negligible and will be used for analyses.") }
    
    if(verbose) { message("Running Wilcoxon tests to test for differences in balances between groups at each non-negligible node.") }
    
    pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(calculated_balances$balances, group1_samples, group2_samples, corr_method=balance_correction, skip_wilcoxon=FALSE)
    
    sig_nodes <- names(calculated_balances$balances)[which(pairwise_node_out$wilcox_corrected_p < balance_p_cutoff)]
    
  } else {
    
    sig_nodes <- significant_nodes
    
    if(any(! sig_nodes %in% phylogeny$node.label)) { stop("Not all sig. nodes are not found in phylogeny.")}
    
    calculated_balances <- list()
    calculated_balances$balances <- tested_balances
    calculated_balances$features <- parallel::mclapply(names(tested_balances),
                                                       lhs_rhs_asvs,
                                                       tree=phylogeny,
                                                       get_node_index=TRUE,
                                                       mc.cores=ncores)
    
    names(calculated_balances$features) <- names(tested_balances)
    
    negligible_nodes_i <- which(! names(phylogeny$node.label) %in% names(tested_balances))
    if(length(negligible_nodes_i) > 0) {
      calculated_balances$negligible_nodes <- phylogeny$node.label[negligible_nodes_i]
    } else {
      calculated_balances$negligible_nodes <- c()
    }
    
    pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(calculated_balances$balances, group1_samples, group2_samples, skip_wilcoxon=TRUE)
  }
  
  if(verbose) { message("Identifying enriched functions at all non-negligible nodes.") }
  
  all_balances_enriched_funcs <- mclapply(names(calculated_balances$balances),
                                          function(x) {
                                            return(node_func_fisher(node = x,
                                                                    in_tree = phylogeny,
                                                                    in_func = func,
                                                                    higher_group=pairwise_node_out$mean_direction[x],
                                                                    pseudocount=1,
                                                                    multiple_test_corr=function_correction))
                                          },
                                          mc.cores = ncores)
  
  names(all_balances_enriched_funcs) <- names(calculated_balances$balances)
  
  if(verbose) { message("Summarizing significant functions across nodes.") }
  
  func_summaries <- summarize_node_enrichment(all_balances_enriched_funcs, sig_nodes, function_p_cutoff)
  
  # Get single DF summarizing the key metrics and print this out.
  func_descrip <- read.table(func_descrip_infile,
                             header=FALSE, sep="\t", row.names=1, stringsAsFactors = FALSE, quote="")
  
  all_func_id <- c()
  for(balance in names(all_balances_enriched_funcs)) {
    all_func_id <- c(all_func_id, rownames(all_balances_enriched_funcs[[balance]]))
  }
  
  all_func_id <- all_func_id[-which(duplicated(all_func_id))]
  
  if(verbose) { message("Creating results dataframe.") }
  summary_df <- data.frame(matrix(NA, nrow=length(all_func_id), ncol=15))
  
  rownames(summary_df) <- all_func_id
  
  colnames(summary_df) <- c("func",
                            "num_sig_nodes_nonenrich",
                            "num_sig_nodes_pos_enrich",
                            "num_sig_nodes_neg_enrich",
                            "num_sig_nodes_not_present",
                            "num_nonsig_nodes_nonenrich",
                            "num_nonsig_nodes_enrich",
                            "num_nonsig_nodes_not_present",
                            "mean_internode_dist_present",
                            "max_internode_dist_present",
                            "mean_internode_dist_neg_enrich",
                            "max_internode_dist_neg_enrich",
                            "mean_internode_dist_pos_enrich",
                            "max_internode_dist_pos_enrich",
                            "description")
  
  rownames(summary_df) <- all_func_id                                  
  summary_df$func <- all_func_id
  summary_df$description <- func_descrip[summary_df$func, "V2"]
  
  phylogeny_node_dists <- dist.nodes(phylogeny)
  
  for(func_id in all_func_id) {
    
    summary_df[func_id, c("num_sig_nodes_nonenrich",
                          "num_sig_nodes_pos_enrich",
                          "num_sig_nodes_neg_enrich",
                          "num_sig_nodes_not_present",
                          "num_nonsig_nodes_nonenrich",
                          "num_nonsig_nodes_enrich",
                          "num_nonsig_nodes_not_present")] <- c(length(func_summaries[[func_id]]$nonenriched_sig_nodes),
                                                                length(func_summaries[[func_id]]$positive_nodes),
                                                                length(func_summaries[[func_id]]$negative_nodes),
                                                                length(func_summaries[[func_id]]$missing_sig_nodes),
                                                                length(func_summaries[[func_id]]$nonenriched_nonsig_nodes),
                                                                length(func_summaries[[func_id]]$enriched_nonsig_nodes),
                                                                length(func_summaries[[func_id]]$missing_nonsig_nodes))
    
    all_nodes_present <- c(func_summaries[[func_id]]$nonenriched_sig_nodes,
                           func_summaries[[func_id]]$positive_nodes,
                           func_summaries[[func_id]]$negative_nodes,
                           func_summaries[[func_id]]$nonenriched_nonsig_nodes,
                           func_summaries[[func_id]]$enriched_nonsig_nodes)
    
    if(max(table(all_nodes_present)) > 1) {
      stop("Node categorized into at least 2 mutually exclusive groups.")
    }
    
    summary_df[func_id, c("mean_internode_dist_present",
                          "max_internode_dist_present",
                          "mean_internode_dist_neg_enrich",
                          "max_internode_dist_neg_enrich",
                          "mean_internode_dist_pos_enrich",
                          "max_internode_dist_pos_enrich")] <- c(internode_mean_max_dist(phy = phylogeny, dist_matrix = phylogeny_node_dists,
                                                                                         node_labels = all_nodes_present),
                                                                 internode_mean_max_dist(phy = phylogeny, dist_matrix = phylogeny_node_dists,
                                                                                         node_labels = func_summaries[[func_id]]$negative_nodes),
                                                                 internode_mean_max_dist(phy = phylogeny, dist_matrix = phylogeny_node_dists,
                                                                                         node_labels = func_summaries[[func_id]]$positive_nodes))
    
  }
  
  sig_nodes_enriched_funcs <- all_balances_enriched_funcs[sig_nodes]
  
  return(list(balances_info=calculated_balances,
              sig_nodes=sig_nodes,
              funcs_per_node=sig_nodes_enriched_funcs,
              df=summary_df,
              out_list=func_summaries,
              tree=phylogeny,
              input_param=input_param))
  
}


parse_taxa_to_df <- function(infile, asvs2keep) {
  
  # Taxonomy files created with commands like this: 
  # tail -n +2 ibd_gevers_2014.otu_table.100.denovo.rdp_assigned | awk '{ print $1 }' > OTU_taxonomy.txt
  disease_taxa <- read.table(file = infile, stringsAsFactors = FALSE, header=FALSE)$V1
  
  # Subset to OTUs in those in input table.
  names(disease_taxa) <- gsub("^.*;d__denovo", "denovo", disease_taxa)
  disease_taxa <- disease_taxa[asvs2keep]
  
  taxa_table <- data.frame(matrix(NA, nrow=length(disease_taxa), ncol=7))
  colnames(taxa_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(taxa_table) <- names(disease_taxa)
  
  taxa_labels <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__", "d__")
  
  for(otu in names(disease_taxa)) {
    
    for(taxon_i in 1:7) {
      str_match <- paste(taxa_labels[taxon_i], ".*;", taxa_labels[taxon_i + 1], sep="")
      otu_taxon <- stringr::str_extract(disease_taxa[otu], str_match)
      otu_taxon <- sub(taxa_labels[taxon_i], "", otu_taxon)
      otu_taxon <- sub(paste(";", taxa_labels[taxon_i + 1], sep=""), "", otu_taxon)
      if(otu_taxon != "") {
        taxa_table[otu, taxon_i] <- otu_taxon
      }
    }
  }
  
  return(taxa_table)
}


determine_min_sibling_group <- function(tree, sibling_tips, ideal_clade_size, min_clade_size) {
  
  # Check that sibling tips are monophyletic.
  if(length(sibling_tips) > 1) {
    ancestral_node <- getMRCA(tree, sibling_tips)
  } else {
    ancestral_node <- which(tree$tip.label == sibling_tips[1])
  }
  
  tips_match_check <- check_node_tips_match_expected(tree = tree, node = ancestral_node, expected_tips = sort(sibling_tips))
  if(! tips_match_check$match) {
    stop(tips_match_check$err)
  }
  
  # Check if sibling tips are already a sufficient group size.
  if(length(sibling_tips) >= ideal_clade_size) {
    return(tree$tip.label[sibling_tips]) 
  }
  
  # If not large enough then loop through all ancestral nodes and return descendants when sufficient tips present.
  for(ancestor in Ancestors(x = tree, node = ancestral_node, type = "all")) {
    sibling_tips <- tree$tip.label[Descendants(x = tree, node = ancestor, type = "tips")[[1]]]
    
    if(length(sibling_tips) >= ideal_clade_size) {
      return(sibling_tips)
    }
  }
  
  # Still return if at least the min set size by this point.
  if(length(sibling_tips) >= min_clade_size) {
    return(sibling_tips)
  }
  
  stop("Error - could not determine an adequately sized sibling group.")
  
}

check_node_tips_match_expected <- function(tree, node, expected_tips, names=TRUE) {
  node_descendants_i <- sort(Descendants(tree, node, type = "tips")[[1]])
  
  if(names) {
    node_descendants <- sort(tree$tip.label[node_descendants_i])
  } else {
    node_descendants <- node_descendants_i
  }
  
  if(! identical(node_descendants, expected_tips)) {
    error_message <- paste("Expected ", length(expected_tips), " tips and found ", length(node_descendants),
                           ", ", length(which(! node_descendants %in% expected_tips)), " of which were not in the expected set.", sep="")
    return(list(match=FALSE, err=error_message))
  }
  
  return(list(match=TRUE))
}

calc_balances <- function(abun_table, lhs_features, rhs_features, pseudocount=NULL) {
  
  if(pseudocount) {
    abun_table <- abun_table + pseudocount
  }

  num_lhs <- length(lhs_features)
  num_rhs <- length(rhs_features)
  
  half1 <- sqrt((num_lhs * num_rhs) / (num_lhs + num_rhs))
  
  sample_col <- colnames(abun_table)
  
  half2 <- c()
  
  for(sample in sample_col) {
    sample_lhs_values <- abun_table[lhs_features, sample]
    sample_rhs_values <- abun_table[rhs_features, sample]

    gmean_lhs <- exp(mean(log(sample_lhs_values)))
    gmean_rhs <- exp(mean(log(sample_rhs_values)))
    
    half2 <- c(half2, log(gmean_lhs / gmean_rhs))
  }
  
  balances_out <- half1 * half2
  
  names(balances_out) <- sample_col
  
  return(balances_out)
}

func_enriched_pathways <- function(funcs, path2funcsets, background_pathway_func_counts, total_funcs) {
  pathway_fisher_tests <- c()
  
  for(pathway in names(path2funcsets)) {
    
    balance_positive_count <- length(which(funcs %in% path2funcsets[[pathway]]))
    
    fisher_test <- fisher.test(matrix(c(balance_positive_count,
                                        length(funcs) - balance_positive_count,
                                        background_pathway_func_counts[pathway],
                                        total_funcs - background_pathway_func_counts[pathway]),
                                      nrow=2))
    
    
    
    pathway_fisher_tests <- c(pathway_fisher_tests, fisher_test$p.value)
    
  }
  pathway_fisher_tests_fdr <- p.adjust(pathway_fisher_tests, "fdr")
  
  names(pathway_fisher_tests_fdr) <- names(path2funcsets)
  
  return(pathway_fisher_tests_fdr)
}

func_pathway_set <- function(pathway2func) {
  
  pathway2funcsets <- list()
  for(pathway in rownames(pathway2func)) {
    
    pathway_funcs <- stringr::str_split(pathway2func[pathway, 1], pattern = ",")[[1]]
    
    if(length(which(pathway_funcs == "")) > 0) {
      pathway_funcs <- pathway_funcs[-which(pathway_funcs == "")]
    }
    
    pathway2funcsets[[pathway]] <- pathway_funcs
  }
  
  return(pathway2funcsets)
  
}

lhs_rhs_asvs <- function(tree, node, get_node_index=FALSE) {

  if(get_node_index) {
    node <- which(tree$node.label == node) + length(tree$tip.label)
  }
  
  node_children <- phangorn::Children(tree, node)

  node_lhs_descendants <- tree$tip.label[phangorn::Descendants(x = tree, node = node_children[1], type = "tips")[[1]]]
  
  node_rhs_descendants <- tree$tip.label[phangorn::Descendants(x = tree, node = node_children[2], type = "tips")[[1]]]
  
  return(list(lhs=node_lhs_descendants,
              rhs=node_rhs_descendants,
              count=c(length(node_lhs_descendants),
                      length(node_rhs_descendants)),
              node_i=node))
}

feature_sets_func_fisher <- function(in_func, feature_set1, feature_set2, pseudocount=NULL, multiple_test_corr="none") {
  
  in_func_set1 <- in_func[feature_set1, ]
  in_func_set2 <- in_func[feature_set2, ]
  
  func_ids_all <- colnames(in_func)
  
  set1_func_missing <- func_ids_all[which(colSums(in_func_set1[, func_ids_all]) == 0)]
  set2_func_missing <- func_ids_all[which(colSums(in_func_set2[, func_ids_all]) == 0)]
  both_func_missing <- set1_func_missing[which(set1_func_missing %in% set2_func_missing)]
  
  if(length(both_func_missing) > 0) {
    in_func_set1 <- in_func_set1[, -which(colnames(in_func_set1) %in% both_func_missing)]
    in_func_set2 <- in_func_set2[, -which(colnames(in_func_set2) %in% both_func_missing)]
  }
  
  func_ids <- colnames(in_func_set1)
  
  fisher_out <- data.frame(matrix(NA, nrow=length(func_ids), ncol=6))
  colnames(fisher_out) <- c("set1_pos", "set1_neg", "set2_pos", "set2_neg", "OR", "P")
  
  rownames(fisher_out) <- func_ids
  
  for(func in func_ids) {
    
    func_count_matrix <- matrix(c(length(which(in_func_set1[, func] > 0)), length(which(in_func_set1[, func] == 0)),
                                  length(which(in_func_set2[, func] > 0)), length(which(in_func_set2[, func] == 0))),
                                nrow=2, ncol=2)
    
    if(pseudocount) {
      func_count_matrix <- func_count_matrix + pseudocount
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
    node_fisher_tests <- feature_sets_func_fisher(in_func = in_func, feature_set1 = node_features$lhs, feature_set2 = node_features$rhs, pseudocount=pseudocount, multiple_test_corr=multiple_test_corr)
  } else if(higher_group == "group2") {
    node_fisher_tests <- feature_sets_func_fisher(in_func = in_func, feature_set1 = node_features$rhs, feature_set2 = node_features$lhs, pseudocount=pseudocount, multiple_test_corr=multiple_test_corr)
  }
  return(node_fisher_tests)
}


calc_OR <- function(exposed_pos, exposed_neg, control_pos, control_neg) {
  return((exposed_pos / exposed_neg) / (control_pos / control_neg))  
}

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



wilcoxon_2group_pvalues <- function(intable, group1_samples, group2_samples) {
  group1_intable <- intable[, group1_samples]
  group2_intable <- intable[, group2_samples]
  
  group1_intable_relab <- data.frame(t(sweep(x = group1_intable, MARGIN = 2, STATS = colSums(group1_intable), FUN = '/')), check.names=FALSE)
  group2_intable_relab <- data.frame(t(sweep(x = group2_intable, MARGIN = 2, STATS = colSums(group2_intable), FUN = '/')), check.names=FALSE)
  
  wilcoxon_p <- c()
  
  for(feature in colnames(group1_intable_relab)) {
    wilcoxon_p <- c(wilcoxon_p, wilcox.test(group1_intable_relab[, feature], group2_intable_relab[, feature])$p.value)
  }
  
  return(wilcoxon_p)
}

compute_tree_node_balances <- function(phylogeny, abun, min_num_tips, ncores=1, pseudocount=1, subset2test=NULL) {
  
  ### Function to perform isomatric log-ratio transformation of feature abundances at each node in the tree.
  ### Will return a list containing the features on the left-hand side (lhs) and right-hand side (rhs) of each
  ### node in a tree ("features") and also the computed balances ("balances").
  
  if(is.null(phylogeny$node.label)) {
   stop("Stopping - input tree does not have any node labels.") 
  }
  
  # Test all nodes unless subset specified.  
  if(! is.null(subset2test)) {
    nodes2test <- subset2test
  } else {
    nodes2test <- phylogeny$node.label
  }

  # Get ASVs on either side of each node.
  node_features <- parallel::mclapply(nodes2test,
                                         lhs_rhs_asvs,
                                         tree=phylogeny,
                                         get_node_index=TRUE,
                                         mc.cores=ncores)
  
  names(node_features) <- nodes2test
  
  negligible_nodes <- sapply(names(node_features),
                                function(x) {
                                  lhs_feat_num <- length(node_features[[x]]$lhs)
                                  rhs_feat_num <- length(node_features[[x]]$rhs)
                                  if((lhs_feat_num < min_num_tips) || (rhs_feat_num < min_num_tips)) {
                                    return(TRUE)
                                  } else {
                                    return(FALSE) 
                                  }
                                })
  
  negligible_nodes_i <- which(negligible_nodes)
  if(length(negligible_nodes_i) > 0) {
    negligible_nodes <- names(negligible_nodes[negligible_nodes_i])
  } else {
    negligible_nodes <- c()
  }
  
  nonnegligible_nodes_i <- which(! names(node_features) %in% negligible_nodes)
  
  if(length(nonnegligible_nodes_i) > 0) {
    nonnegligible_nodes <- names(node_features)[nonnegligible_nodes_i]
  } else {
    stop("Stopping - no non-negligible nodes remain after filtering based on mininum number of tips of left and right-hand side of each node.")
  }
  
  # Calculate balances at each node.
  balance_calc <- mclapply(nonnegligible_nodes,
                           function(x) {
                             return(calc_balances(abun_table=abun,
                                                  lhs_features=node_features[[x]]$lhs,
                                                  rhs_features=node_features[[x]]$rhs,
                                                  pseudocount=pseudocount))
                           },
                           mc.cores=ncores)
  
  names(balance_calc) <- nonnegligible_nodes
  
  return(list(features=node_features, balances=balance_calc, negligible_nodes=negligible_nodes))
  
}

balance_taxa_name <- function(lhs_features, rhs_features, taxa, threshold=0.75) {
  
  if(threshold <= 0.5 | threshold > 1) {
    stop("The set threshold needs to be > 0.5 and <= 1.")
  }
  
  taxon_lhs <- feature_consensus_taxon(taxa, lhs_features, threshold)
  taxon_rhs <- feature_consensus_taxon(taxa, rhs_features, threshold)
  
  return(paste(taxon_lhs, taxon_rhs, sep="/"))
}

feature_consensus_taxon <- function(taxa, features, threshold) {
  
  taxa_subset <- taxa[features, ]
  
  for(i in ncol(taxa_subset):1) {
    
    taxa_subset_table <- table(taxa_subset[, i, drop=TRUE])
    
    if(length(taxa_subset_table) == 0) { next }
    
    if((max(taxa_subset_table) / length(features)) >= threshold) {
      return(paste(colnames(taxa_subset)[i], names(taxa_subset_table)[which(taxa_subset_table == max(taxa_subset_table))], sep="_"))
    }
  }
  
  return("Unclear")
}

summarize_node_enrichment <- function(enriched_funcs, sig_nodes, func_p_cutoff) {
  # For each function that is significant at least once get:
  #   - Names of significant nodes where the function is present.
  #   - Names of positively and negatively enriched significant nodes
  #   - Names of nonsignificant nodes where the function would have been significant.
  #   - Names of nonsignificant nodes where the function was present.
  
  # First get all unique functions.
  all_func_id <- c()
  for(node in names(enriched_funcs)) {
    all_func_id <- c(all_func_id, rownames(enriched_funcs[[node]]))
  }
  
  all_func_id <- all_func_id[-which(duplicated(all_func_id))]
  
  func_summaries <- list()
  
  # Loop through each function and get breakdown of contributing node names.
  for(func_id in all_func_id) {
    
    func_summaries[[func_id]]$positive_nodes <- c()
    func_summaries[[func_id]]$negative_nodes <- c()
    func_summaries[[func_id]]$nonenriched_sig_nodes <- c()
    func_summaries[[func_id]]$nonenriched_nonsig_nodes <- c()
    func_summaries[[func_id]]$enriched_nonsig_nodes <- c()
    func_summaries[[func_id]]$missing_sig_nodes <- c()
    func_summaries[[func_id]]$missing_nonsig_nodes <- c()
    
    for(node in names(enriched_funcs)) {
      if(func_id %in% rownames(enriched_funcs[[node]])) {
        
        if(node %in% sig_nodes) {
          # If this node was significant then categorize it as either positive, negative, or nonenriched.
          if(as.numeric(enriched_funcs[[node]][func_id, "P_corr"]) < func_p_cutoff) {
            
            if(as.numeric(enriched_funcs[[node]][func_id, "OR"]) > 1) {
              func_summaries[[func_id]]$positive_nodes <- c(func_summaries[[func_id]]$positive_nodes, node)
            } else if(as.numeric(enriched_funcs[[node]][func_id, "OR"]) < 1) {
              func_summaries[[func_id]]$negative_nodes <- c(func_summaries[[func_id]]$negative_nodes, node)
            } else {
              print(enriched_funcs[[node]][func_id, ])
              stop("Significant function but OR for above info not different from 1?!")
            }
          } else {
            func_summaries[[func_id]]$nonenriched_sig_nodes <- c(func_summaries[[func_id]]$nonenriched_sig_nodes, node)
          }
        } else {
          # Since this node was NOT significant then categorize it as either enriched or nonenriched.
          if(as.numeric(enriched_funcs[[node]][func_id, "P_corr"]) < func_p_cutoff) {
            func_summaries[[func_id]]$enriched_nonsig_nodes <- c(func_summaries[[func_id]]$enriched_nonsig_nodes, node)
          } else {
            func_summaries[[func_id]]$nonenriched_nonsig_nodes <- c(func_summaries[[func_id]]$nonenriched_nonsig_nodes, node)
          }
        }
      } else {
        if(node %in% sig_nodes) {
          func_summaries[[func_id]]$missing_sig_nodes <- c(func_summaries[[func_id]]$missing_sig_nodes, node)
        } else {
          func_summaries[[func_id]]$missing_nonsig_nodes <- c(func_summaries[[func_id]]$missing_nonsig_nodes, node)
        }
      }
    }
  }
  return(func_summaries)
}

internode_mean_max_dist <- function(phy, dist_matrix, node_labels) {
  
  # If fewer than 2 nodes in set then return NA.
  if(length(node_labels) <= 1) {
    missing_val <- c(NA, NA)
    names(missing_val) <- c("mean", "max")
    return(missing_val)
  }
  
  nodes_i <- which(phy$node.label %in% node_labels) + length(phy$tip.label)
  
  if(length(nodes_i) != length(node_labels)) {
    stop(paste("Wrong number of node labels match input set. Looking for ",
               paste(node_labels, collapse = ", "),
               " and found indices ",
               paste(as.character(nodes_i), collapse=", "),
               ".",
               sep=""))
  }
  
  dist_matrix_subset <- as.dist(dist_matrix[nodes_i, nodes_i])
  
  internode_dist_metrics <- c(mean(dist_matrix_subset), max(dist_matrix_subset))
  names(internode_dist_metrics) <- c("mean", "max")
  
  return(internode_dist_metrics)
  
}

prep_tree_sig_balances <- function(in_list, taxa_table) {
  
  all_sig_balances <- c()
  for(func in names(in_list$out_list$balances)) {
    all_sig_balances <- c(all_sig_balances,
                          in_list$out_list$balances[[func]]$positive_balances,
                          in_list$out_list$balances[[func]]$negative_balances)
  }
  
  all_sig_balances <- all_sig_balances[-which(duplicated(all_sig_balances))]
  
  sig_node_taxa <- list()
  
  for(node in all_sig_balances) {
    sig_node_taxa[[node]] <- balance_taxa_name(lhs_features = in_list$balances_info$features[[node]]$lhs,
                                               rhs_features = in_list$balances_info$features[[node]]$rhs,
                                               taxa = taxa_table)
  }
  
  tree_sig_subset <- in_list$tree
  nodes2ignore <- which(! tree_sig_subset$node.label %in% all_sig_balances)
  for(node in all_sig_balances) {
    tree_sig_subset$node.label[which(tree_sig_subset$node.label == node)] <- sig_node_taxa[[node]]
  }
  tree_sig_subset$node.label[nodes2ignore] <- ""
  
  all_sig_balances_i <- which(tree_sig_subset$node.label != "") + length(tree_sig_subset$tip.label)
  
  return(list(prepped_tree=tree_sig_subset, nodes2plot=all_sig_balances_i))
}

prep_tree_balances_func <- function(in_list, focal_func, taxa_table) {
  
  all_sig_balances <- c()
  for(func in names(in_list$out_list$balances)) {
    all_sig_balances <- c(all_sig_balances,
                          in_list$out_list$balances[[func]]$positive_balances,
                          in_list$out_list$balances[[func]]$negative_balances)
  }
  
  all_sig_balances <- all_sig_balances[-which(duplicated(all_sig_balances))]
  
  sig_node_taxa <- list()
  
  for(node in all_sig_balances) {
    sig_node_taxa[[node]] <- balance_taxa_name(lhs_features = in_list$balances_info$features[[node]]$lhs,
                                               rhs_features = in_list$balances_info$features[[node]]$rhs,
                                               taxa = taxa_table)
  }
  
  tree_sig_subset <- in_list$tree
  nodes2ignore <- which(! tree_sig_subset$node.label %in% all_sig_balances)
  for(node in all_sig_balances) {
    tree_sig_subset$node.label[which(tree_sig_subset$node.label == node)] <- sig_node_taxa[[node]]
  }
  tree_sig_subset$node.label[nodes2ignore] <- ""
  
  enriched_balances <- c()
  depleted_balances <- c()
  
  for(balance in all_sig_balances) {
    
    if(! focal_func %in% rownames(in_list$funcs_per_balance[[balance]])) {
      next
    }
    
    if(in_list$funcs_per_balance[[balance]][focal_func, "P"] < 0.05) {
      if(in_list$funcs_per_balance[[balance]][focal_func, "OR"] > 1) {
        enriched_balances <- c(enriched_balances, balance)
      } else if(in_list$funcs_per_balance[[balance]][focal_func, "OR"] < 1) {
        depleted_balances <- c(depleted_balances, balance)
      }
    }
  }
  
  enriched_balances_i <- which(in_list$tree$node.label %in% enriched_balances) + length(in_list$tree$tip.label)
  
  depleted_balances_i <- which(in_list$tree$node.label %in% depleted_balances) + length(in_list$tree$tip.label)
  
  all_sig_balances_minus_enriched_depleted <- all_sig_balances[which(! all_sig_balances %in% enriched_balances)]
  all_sig_balances_minus_enriched_depleted <- all_sig_balances_minus_enriched_depleted[which(! all_sig_balances_minus_enriched_depleted %in% depleted_balances)]
  
  nonsig_balances_i <- which(in_list$tree$node.label %in% all_sig_balances_minus_enriched_depleted) + length(in_list$tree$tip.label)
  
  return(list(prepped_tree=tree_sig_subset, enriched=enriched_balances_i, depleted=depleted_balances_i, nonsig=nonsig_balances_i))
}

categorize_gene_families <- function(in_list, min_sig_balances=4, min_prop_enriched_balances=0.3, max_opposite_enriched=0.1) {
  
  in_list$df$prop_sig_balances_pos_enrich <- in_list$df$num_sig_balances_pos_enrich / in_list$df$num_sig_balances_present
  
  in_list$df$prop_sig_balances_neg_enrich <- in_list$df$num_sig_balances_neg_enrich / in_list$df$num_sig_balances_present
  
  
  in_list_pos_func <- in_list$df[which((in_list$df$num_sig_balances_present >= min_sig_balances) & (in_list$df$prop_sig_balances_pos_enrich >= min_prop_enriched_balances) & (in_list$df$prop_sig_balances_neg_enrich <max_opposite_enriched)), ]
  
  in_list_neg_func <- in_list$df[which((in_list$df$num_sig_balances_present >= min_sig_balances) & (in_list$df$prop_sig_balances_neg_enrich >= min_prop_enriched_balances) & (in_list$df$prop_sig_balances_pos_enrich < max_opposite_enriched)), ]
  
  in_list_mixed_func <- in_list$df[which((in_list$df$num_sig_balances_present >= min_sig_balances) & (in_list$df$prop_sig_balances_neg_enrich > max_opposite_enriched) & (in_list$df$prop_sig_balances_pos_enrich > max_opposite_enriched)), ]
  
  in_list_nonsig_func <- in_list$df[which((in_list$df$num_sig_balances_present > 0) & (in_list$df$prop_sig_balances_neg_enrich == 0) & (in_list$df$prop_sig_balances_pos_enrich == 0)), ]
  
  in_list_other_func <- in_list$df[which(! in_list$df$func %in% c(in_list_pos_func$func, in_list_neg_func$func, in_list_mixed_func$func, in_list_nonsig_func$func)), ]
  
  return(list(pos=in_list_pos_func, neg=in_list_neg_func, mixed=in_list_mixed_func, nonsig=in_list_nonsig_func, other=in_list_other_func))
}

parse_sig_pathway_names <- function(pathway_list) {
  
  out_names <- c()
  
  for(pathway in names(pathway_list)) {
    out_names <- c(out_names, paste(pathway, pathway_list[[pathway]]$descrip, sep=" - "))
  }
  return(out_names)
}


pairwise_mean_direction_and_wilcoxon <- function(in_list, group1, group2, corr_method="BH", ncores=1, skip_wilcoxon=FALSE) {
  
  wilcox_corrected_p <- NULL
  wilcox_raw_p <- NULL
  wilcox_output <- list()
  mean_direction <- c()
  
  if(skip_wilcoxon) {
    result <- parallel::mclapply(names(in_list), function(x) {
      
      group1_mean <- mean(in_list[[x]][group1])
      group2_mean <- mean(in_list[[x]][group2])
      
      if(group1_mean > group2_mean) {
        mean_direction <- c(mean_direction, "group1")
      } else if(group1_mean < group2_mean) {
        mean_direction <- c(mean_direction, "group2")
      } else if(group1_mean == group2_mean) {
        mean_direction <- c(mean_direction, "same")
        warning(paste("The calculated means are exactly the same for each group for test ", x, ", which likely indicates a problem.", sep=""))
      }
      
      return(mean_direction)
    }, mc.cores=ncores)
    
    for(i in 1:length(result)) {
      mean_direction <- c(mean_direction, result[[i]])
    }
    
    
  } else {
    result <- parallel::mclapply(names(in_list), function(x) {
      wilcox_out <- wilcox.test(in_list[[x]][group1], in_list[[x]][group2], exact=FALSE)
      
      group1_mean <- mean(in_list[[x]][group1])
      group2_mean <- mean(in_list[[x]][group2])
      
      if(group1_mean > group2_mean) {
        mean_direction <- c(mean_direction, "group1")
      } else if(group1_mean < group2_mean) {
        mean_direction <- c(mean_direction, "group2")
      } else if(group1_mean == group2_mean) {
        mean_direction <- c(mean_direction, "same")
        warning(paste("The calculated means are exactly the same for each group for test ", x, ", which likely indicates a problem.", sep=""))
      }
      
      return(list(wilcox_out=wilcox_out, mean_direction=mean_direction))
    }, mc.cores=ncores)
    
    wilcox_raw_p <- c()
    
    for(i in 1:length(result)) {
      wilcox_raw_p <- c(wilcox_raw_p, result[[i]]$wilcox_out$p.value)
      wilcox_output[[i]] <- result[[i]]$wilcox_out
      mean_direction <- c(mean_direction, result[[i]]$mean_direction)
    }
    
    wilcox_corrected_p <- p.adjust(p = wilcox_raw_p, method = corr_method)
    
    names(wilcox_raw_p) <- names(in_list)
    names(wilcox_corrected_p) <- names(in_list)
    names(wilcox_output) <- names(in_list)
  }
  
  names(mean_direction) <- names(in_list)
  
  return(list(mean_direction=mean_direction, wilcox_raw_p=wilcox_raw_p, wilcox_corrected_p=wilcox_corrected_p, wilcox_output=wilcox_output))
}


filter_rare_table_cols <- function(in_tab, min_nonzero_count, min_nonzero_prop, verbose=TRUE) {
  nonzero_counts <- colSums(in_tab > 0)
  col2remove <- which((nonzero_counts < min_nonzero_count) | (nonzero_counts / ncol(in_tab) < min_nonzero_prop))
  if(length(col2remove) > 0) {
    if(verbose) { message(paste("Filtering", as.character(length(col2remove)), "rare functions from input function table.")) }
    in_tab <- in_tab[, -col2remove]
  }
  return(in_tab)
}


check_two_group_balance_args <- function(abun,
                                         func,
                                         phylogeny,
                                         taxa,
                                         group1_samples,
                                         group2_samples,
                                         ncores,
                                         pseudocount,
                                         significant_nodes,
                                         tested_balances,
                                         min_num_tips,
                                         min_func_instances,
                                         min_func_prop,
                                         name_threshold,
                                         balance_p_cutoff,
                                         balance_correction,
                                         function_p_cutoff,
                                         function_correction,
                                         func_descrip_infile,
                                         verbose) {
  
  if(((is.null(significant_nodes)) && (! is.null(tested_balances))) || ((! is.null(significant_nodes)) && (is.null(tested_balances)))) {
    stop("Stopping - arguments significant_nodes and tested_balances either both need to be given or neither should be specified.")
  }
  
  if((! is.null(significant_nodes)) && (length(significant_nodes) == 0)) { stop("Stopping - vector specified for significant_nodes argument is empty.") }
  
  if(class(abun) != "data.frame") { stop("Stopping - argument abun needs to be of the class data.frame.") }
  if(class(func) != "data.frame") { stop("Stopping - argument func needs to be of the class data.frame.") }
  if(class(taxa) != "data.frame") { stop("Stopping - argument taxa needs to be of the class data.frame.") }
  if(class(phylogeny) != "phylo") { stop("Stopping - argument phylo needs to be of the class phylo.") }
  
  if(class(group1_samples) != "character") { stop("Stopping - argument group1_samples needs to be of the class character.") }
  if(class(group2_samples) != "character") { stop("Stopping - argument group2_samples needs to be of the class character.") }
  
  if(length(group1_samples) == 0) { stop("Stopping - argument group1_samples is of length 0.") }
  if(length(group2_samples) == 0) { stop("Stopping - argument group2_samples is of length 0.") }
  
  if(length(which(group1_samples %in% colnames(abun))) != length(group1_samples)) { stop("Stopping - not all group1_samples match columns of abun argument.") }
  if(length(which(group2_samples %in% colnames(abun))) != length(group2_samples)) { stop("Stopping - not all group2_samples match columns of abun argument.") }
  
  if((class(ncores) != "integer") && (class(ncores) != "numeric")) { stop("Stopping - ncores argument needs to be of class numeric or integer.") }
  if(ncores <= 0) { stop("Stopping - ncores argument needs to be higher than 0.") }
  
  if((class(pseudocount) != "integer") && (class(pseudocount) != "numeric")) { stop("Stopping - pseudocount argument needs to be of class numeric or integer.") }
  if(pseudocount < 0) { stop("Stopping - pseudocount argument cannot be lower than 0.") }
  
  if(min_num_tips > length(phylogeny$tip.label) / 2) { stop("Stopping - the min_num_tips argument cannot be higher than half of the total number of tips.") }
  
  if((min_func_prop < 0) || (min_func_prop > 1)) { stop("Stopping - the min_func_prop argument must be between 0 and 1.") }
  if((balance_p_cutoff < 0) || (balance_p_cutoff > 1)) { stop("Stopping - the balance_p_cutoff argument must be between 0 and 1.") }
  if((function_p_cutoff < 0) || (function_p_cutoff > 1)) { stop("Stopping - the function_p_cutoff argument must be between 0 and 1.") }
  if((name_threshold < 0) || (name_threshold > 1)) { stop("Stopping - the name_threshold argument must be between 0 and 1.") }
  
  if((min_func_instances < 0) || (min_func_instances > ncol(func))) { stop("Stopping - the min_func_instances argument must be between 0 and 1.") }
  
  if(! function_correction %in% p.adjust.methods) { stop("Stopping - function_correction argument needs to be found in p.adjust.methods.") }
  if(! balance_correction %in% p.adjust.methods) { stop("Stopping - balance_correction argument needs to be found in p.adjust.methods.") }
  
  if(! file.exists(func_descrip_infile)) { stop("Stopping - file corresponding to func_descrip_infile argument not found.") } 
  
  if(! is.logical(verbose)) { stop("Stopping - verbose argument needs to be TRUE or FALSE.") }
  
  return(list(group1_samples=group1_samples, group2_samples=group2_samples,
              ncores=ncores, pseudocount=pseudocount, significant_nodes=significant_nodes,
              min_num_tips=min_num_tips, min_func_instances=min_func_instances,
              min_func_prop=min_func_prop, name_threshold=name_threshold,
              balance_p_cutoff=balance_p_cutoff, balance_correction=balance_correction,
              function_p_cutoff=function_p_cutoff, function_correction=function_correction,
              func_descrip_infile=func_descrip_infile, verbose=verbose))
}

outlier_func_enriched_pathways <- function(summary_df,
                                           pos_cutoff,
                                           neg_cutoff,
                                           p_corr_method="BH",
                                           corr_P_cutoff=0.05,
                                           pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                           pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv") {
  
  
  if(! file.exists(pathway2func_map)) { stop("Stopping - file corresponding to pathway2func_map argument not found.") } 
  if(! file.exists(pathway_descrip_infile)) { stop("Stopping - file corresponding to pathway_descrip_infile argument not found.") } 
  
  # Prep pathway mapfiles.
  pathway_descrip <- read.table(pathway_descrip_infile,
                                header=FALSE, sep="\t", row.names=1, quote="", comment.char="", stringsAsFactors = FALSE)
  
  pathway2func <- read.table(file = pathway2func_map,
                             sep="\t", header=FALSE, stringsAsFactors = FALSE, row.names=1)
  pathway2funcsets <- list()
  
  
  for(pathway in rownames(pathway2func)) {
    
    pathway_funcs <- stringr::str_split(pathway2func[pathway, 1], pattern = ",")[[1]]
    
    if(length(which(pathway_funcs == "")) > 0) {
      pathway_funcs <- pathway_funcs[-which(pathway_funcs == "")]
    }
    
    if(length(which(! pathway_funcs %in% summary_df$func)) > 0) {
      pathway_funcs <- pathway_funcs[-which(! pathway_funcs %in% summary_df$func)]
    }
    
    pathway2funcsets[[pathway]] <- pathway_funcs
  }
  
  # Get background of how many functions in this dataset are in each pathway.
  total_funcs <- length(summary_df$func)
  background_pathway_func_counts <- c()
  
  for(pathway in names(pathway2funcsets)) {
    func_set <- pathway2funcsets[[pathway]]
    func_set <- func_set[which(func_set %in% summary_df$func)]
    
    if(length(func_set) <= 2) {
      background_pathway_func_counts <- c(background_pathway_func_counts, NA)
    } else {
      background_pathway_func_counts <- c(background_pathway_func_counts, length(func_set))
    }
  }
  
  names(background_pathway_func_counts) <- names(pathway2funcsets)
  
  if(length(which(is.na(background_pathway_func_counts))) > 0) {
    background_pathway_func_counts <- background_pathway_func_counts[-which(is.na(background_pathway_func_counts))]
    pathway2funcsets <- pathway2funcsets[names(background_pathway_func_counts)]
  }
  
  func_subsets_to_test <- list(up=summary_df$func[which(summary_df$num_sig_nodes_pos_enrich >= pos_cutoff)],
                               down=summary_df$func[which(summary_df$num_sig_nodes_neg_enrich >= neg_cutoff)])
  
  enriched_pathways <- list()
  
  for(n in names(func_subsets_to_test)) {
    pathway_fisher_tests <- c()
    pathway_ORs <- c()
    pathway_raw_counts <- list()
    enriched_pathways[[n]] <- list()
    
    for(pathway in names(pathway2funcsets)) {
      positive_count <- length(which(func_subsets_to_test[[n]] %in% pathway2funcsets[[pathway]]))
      
      pathway_raw_counts[[pathway]] <- matrix(c(positive_count,
                                                length(func_subsets_to_test[[n]]) - positive_count,
                                                background_pathway_func_counts[pathway],
                                                total_funcs - background_pathway_func_counts[pathway]),
                                              nrow=2)
      
      fisher_test <- fisher.test(pathway_raw_counts[[pathway]], alternative="greater")
      
      pathway_fisher_tests <- c(pathway_fisher_tests, fisher_test$p.value)
      
      pathway_ORs <- c(pathway_ORs, calc_OR(exposed_pos = positive_count,
                                            exposed_neg = length(func_subsets_to_test[[n]]) - positive_count,
                                            control_pos = background_pathway_func_counts[pathway],
                                            control_neg = total_funcs - background_pathway_func_counts[pathway]))
      
    }
    
    pathway_fisher_tests_corr <- p.adjust(pathway_fisher_tests, p_corr_method)
    
    for(i in which(pathway_fisher_tests_corr < corr_P_cutoff)) {
      pathway_id <- names(pathway2funcsets)[i]
      enriched_pathways[[n]][[pathway_id]] <- list(P=pathway_fisher_tests[i],
                                                   P_corr=pathway_fisher_tests_corr[i],
                                                   OR=pathway_ORs[i],
                                                   raw_counts=pathway_raw_counts[[pathway_id]],
                                                   descrip=pathway_descrip[pathway_id, "V2"])
    }
  }
  
  return(enriched_pathways)
}