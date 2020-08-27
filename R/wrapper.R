#' Two-group POMS pipeline
#'
#' This function will identify significant nodes based on Wilcoxon tests
#' and then identify enriched functions.
#'
#' @param abun Dataframe of MAG or ASV abundances
#' @return A matrix of the infile
#' @export
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

#' @export
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

