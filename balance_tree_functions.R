library(ape)
library(ALDEx2)
library(corncob)
library(fastcluster)
library(grid)
library(parallel)
library(phangorn)
library(stringr)
library(VennDiagram)

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
      otu_taxon <- str_extract(disease_taxa[otu], str_match)
      otu_taxon <- sub(taxa_labels[taxon_i], "", otu_taxon)
      otu_taxon <- sub(paste(";", taxa_labels[taxon_i + 1], sep=""), "", otu_taxon)
      if(otu_taxon != "") {
        taxa_table[otu, taxon_i] <- otu_taxon
      }
    }
  }
  
  return(taxa_table)
}

determine_balance_groups <- function(func, func_by_asv, asv_dist, ideal_sibling_size, min_sibling_size) {
  
  func_pos_asvs <- which(func_by_asv[, func] > 0)
  
  focal_asvs <- rownames(func_by_asv)[func_pos_asvs]
  
  asv_dist[func_pos_asvs, func_pos_asvs] <- 0
  
  colnames(asv_dist) <- rownames(func_by_asv)
  rownames(asv_dist) <- rownames(func_by_asv)
  
  asv_dist_phylo <- as.phylo(fastcluster::hclust(as.dist(asv_dist)))
  asv_dist_phylo <- makeNodeLabel(asv_dist_phylo, method="number", prefix='n')
  
  if(length(focal_asvs) > 1) {
    func_grouping_mrca <- getMRCA(asv_dist_phylo, focal_asvs)
  } else {
    func_grouping_mrca <- which(asv_dist_phylo$tip.label == focal_asvs[1])
  }
  
  # Double-check that only descendants of MRCA are the tips contributing the same function.
  tips_match_check <- check_node_tips_match_expected(tree = asv_dist_phylo, node = func_grouping_mrca, expected_tips = sort(focal_asvs))
  if(! tips_match_check$match) {
    stop(tips_match_check$err)
  }
  
  func_grouping_mrca_sibling <- Siblings(x = asv_dist_phylo, node = func_grouping_mrca)
  
  sibling_asvs <- asv_dist_phylo$tip.label[Descendants(asv_dist_phylo, func_grouping_mrca_sibling, type = "tips")[[1]]]
  
  if(length(sibling_asvs) < ideal_sibling_size) {
    asv_dist_phylo <- drop.tip(phy = asv_dist_phylo, tip = func_pos_asvs, trim.internal = TRUE)
    sibling_asvs <- determine_min_sibling_group(tree=asv_dist_phylo, sibling_tips=sibling_asvs, ideal_clade_size=ideal_sibling_size, min_clade_size=min_sibling_size)
  }
  
  return(list(focal=focal_asvs, sibling=sibling_asvs))
  
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
    
    pathway_funcs <- str_split(pathway2func[pathway, 1], pattern = ",")[[1]]
    
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

threeWayVennWrapper <- function(set1, set2, set3,
                                labels=c("cat1", "cat2", "cat3"),
                                colours=c("#009E73", "#E69F00", "#56B4E9")) {
  
  set1_count <- length(set1)
  set2_count <- length(set2)
  set3_count <- length(set3)
  set1_2_count <- length(which(set1 %in% set2))
  set2_3_count <- length(which(set2 %in% set3))
  set1_3_count <- length(which(set1 %in% set3))
  set1_2_3_count_TMP <- set1[which(set1 %in% set2)]
  set1_2_3_count <- length(set1_2_3_count_TMP[which(set1_2_3_count_TMP %in% set3)])
  
  grid.newpage()
  
  venn_out <- draw.triple.venn(area1=set1_count,
                               area2=set2_count,
                               area3=set3_count,
                               n12 = set1_2_count,
                               n23 = set2_3_count,
                               n13 = set1_3_count,
                               n123 = set1_2_3_count,
                               category = labels,
                               scaled=TRUE,
                               fill = colours,
                               cex=rep(2, 7),
                               cat.cex=rep(2, 3))
  
  return(venn_out)
  
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

identify_enriched_pathways <- function(func, func_subsets_to_test, pathway2func_map, pathway_descrip_infile) {
  
  # Prep pathway mapfiles.
  pathway_descrip <- read.table(pathway_descrip_infile,
                                header=FALSE, sep="\t", row.names=1, quote="", comment.char="", stringsAsFactors = FALSE)
  
  pathway2func <- read.table(file = pathway2func_map,
                             sep="\t", header=FALSE, stringsAsFactors = FALSE, row.names=1)
  pathway2funcsets <- list()
  for(pathway in rownames(pathway2func)) {
    
    pathway_funcs <- str_split(pathway2func[pathway, 1], pattern = ",")[[1]]
    
    if(length(which(pathway_funcs == "")) > 0) {
      pathway_funcs <- pathway_funcs[-which(pathway_funcs == "")]
    }
    
    pathway2funcsets[[pathway]] <- pathway_funcs
  }
  
  # Get background of how many functions across ASVs in this dataset are in each pathway.
  total_funcs <- sum(colSums(func))
  background_pathway_func_counts <- c()
  
  for(pathway in names(pathway2funcsets)) {
    func_set <- pathway2funcsets[[pathway]]
    func_set <- func_set[which(func_set %in% colnames(func))]
    
    if(length(func_set) <= 2) {
      background_pathway_func_counts <- c(background_pathway_func_counts, NA) 
    } else {
      background_pathway_func_counts <- c(background_pathway_func_counts, sum(func[, func_set]))
    }
  }
  
  names(background_pathway_func_counts) <- names(pathway2funcsets)
  background_pathway_func_counts <- background_pathway_func_counts[-which(is.na(background_pathway_func_counts))]
  pathway2funcsets <- pathway2funcsets[names(background_pathway_func_counts)]
  
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
      
      fisher_test <- fisher.test(matrix(c(positive_count,
                                          length(func_subsets_to_test[[n]]) - positive_count,
                                          background_pathway_func_counts[pathway],
                                          total_funcs - background_pathway_func_counts[pathway]),
                                        nrow=2), alternative="greater")
      
      pathway_fisher_tests <- c(pathway_fisher_tests, fisher_test$p.value)
      
      pathway_ORs <- c(pathway_ORs, calc_OR(exposed_pos = positive_count,
                                            exposed_neg = length(func_subsets_to_test[[n]]) - positive_count,
                                            control_pos = background_pathway_func_counts[pathway],
                                            control_neg = total_funcs - background_pathway_func_counts[pathway]))
      
    }
    
    pathway_fisher_tests_BY <- p.adjust(pathway_fisher_tests, "BY")
    
    for(i in which(pathway_fisher_tests_BY < 0.05)) {
      pathway_id <- names(pathway2funcsets)[i]
      enriched_pathways[[n]][[pathway_id]] <- list(P=pathway_fisher_tests[i],
                                                   BY=pathway_fisher_tests_BY[i],
                                                   OR=pathway_ORs[i],
                                                   raw_counts=pathway_raw_counts[[pathway_id]],
                                                   descrip=pathway_descrip[pathway_id, "V2"])
    }
  }
  
  return(enriched_pathways)
}

calc_func_abun <- function(in_abun, in_func, ncores=1) {
  
  out_df <- data.frame(matrix(NA, nrow=ncol(in_func), ncol=ncol(in_abun)))
  colnames(out_df) <- colnames(in_abun)
  rownames(out_df) <- colnames(in_func)
  
  out_sample_func_abun <- mclapply(colnames(in_abun), function(x) { return(colSums(in_abun[, x] * in_func)) }, mc.cores=ncores)
  names(out_sample_func_abun) <- colnames(in_abun)
  
  for(sample in colnames(in_abun)) {
    out_df[, sample] <- out_sample_func_abun[[sample]]
  }
  
  return(out_df)
  
}

# Run other differential abundance tools.
run_2group_ALDEx2 <- function(in_table, group1_samples, group2_samples) {
  in_table <- round(in_table[, c(group1_samples, group2_samples)])
  return(aldex(reads = in_table, conditions=c(rep("group1", length(group1_samples)), rep("group2", length(group2_samples)))))
}

run_2group_corncob <- function(in_table, group1_samples, group2_samples) {
  
  in_table <- in_table[, c(group1_samples, group2_samples)]
  
  sample_meta <- data.frame(grouping=c(rep("group1", length(group1_samples)),
                                       rep("group2", length(group2_samples))))
  
  rownames(sample_meta) <- colnames(in_table)
  
  return(differentialTest(formula = . ~ grouping,
                          phi.formula = ~ grouping,
                          formula_null = ~ 1,
                          phi.formula_null = ~ grouping,
                          data = in_table,
                          test = "Wald",
                          boot = FALSE,
                          sample_data = sample_meta,
                          taxa_are_rows = TRUE))
}

corncob_determine_sig_sets <- function(in_table, corncob_out, group1_samples, group2_samples) {
  
  corncob_out_all_sig <- names(corncob_out$BY)[which(corncob_out$BY < 0.05)]
  
  corncob_out_group1_higher <- c()
  corncob_out_group1_lower <- c()
  
  for(func in corncob_out_all_sig) {
    mean_diff <- mean(as.numeric(in_table[func, group1_samples])) - mean(as.numeric(in_table[func, group2_samples]))
    if(mean_diff > 0) {
      corncob_out_group1_higher <- c(corncob_out_group1_higher, func)
    } else if(mean_diff < 0) {
      corncob_out_group1_lower <- c(corncob_out_group1_lower, func)
    }
  }
  
  return(list(all_sig=corncob_out_all_sig, group1_higher=corncob_out_group1_higher, group1_lower=corncob_out_group1_lower))
  
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

compute_tree_node_balances <- function(phylogeny, abun, ncores=1, pseudocount=1) {
  
  ### Function to perform isomatric log-ratio transformation of feature abundances at each node in the tree.
  ### Will return a list containing the features on the left-hand side (lhs) and right-hand side (rhs) of each
  ### node in a tree ("features") and also the computed balances ("balances").
  
  # Get ASVs on either side of each node.
  balance_features <- parallel::mclapply(phylogeny$node.label,
                                         lhs_rhs_asvs,
                                         tree=phylogeny,
                                         get_node_index=TRUE,
                                         mc.cores=ncores)
  
  names(balance_features) <- phylogeny$node.label
  
  # Calculate balances at each node.
  balance_calc <- mclapply(names(balance_features),
                           function(x) {
                             return(calc_balances(abun_table=abun,
                                                  lhs_features=balance_features[[x]]$lhs,
                                                  rhs_features=balance_features[[x]]$rhs,
                                                  pseudocount=pseudocount))
                           },
                           mc.cores=ncores)
  
  names(balance_calc) <- names(balance_features)
  
  return(list(features=balance_features, balances=balance_calc))
  
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

summarize_balance_enrichment <- function(enriched_funcs, sig_balances, func_p_cutoff) {
  # For each function that is significant at least once get:
  #   - Names of significant balances where the function is present.
  #   - Names of positively and negatively enriched significant balances.
  #   - Names of nonsignificant balances where the function would have been significant.
  #   - Names of nonsignificant balances where the function was present.
  
  # First get all unique functions.
  all_func_id <- c()
  for(balance in names(enriched_funcs)) {
    all_func_id <- c(all_func_id, rownames(enriched_funcs[[balance]]))
  }
  
  all_func_id <- all_func_id[-which(duplicated(all_func_id))]
  
  func_summaries <- list()
  func_summaries$balances <- list()
  
  # Loop through each function and get breakdown of contributing balance names.
  for(func_id in all_func_id) {
    
    func_summaries$balances[[func_id]]$positive_balances <- c()
    func_summaries$balances[[func_id]]$negative_balances <- c()
    func_summaries$balances[[func_id]]$nonenriched_sig_balances <- c()
    func_summaries$balances[[func_id]]$nonenriched_nonsig_balances <- c()
    func_summaries$balances[[func_id]]$enriched_nonsig_balances <- c()
    func_summaries$balances[[func_id]]$missing_sig_balances <- c()
    func_summaries$balances[[func_id]]$missing_nonsig_balances <- c()
    
    for(balance in names(enriched_funcs)) {
      if(func_id %in% rownames(enriched_funcs[[balance]])) {
        
        if(balance %in% sig_balances) {
          # If this balance was significant then categorize it as either positive, negative, or nonenriched.
          if(as.numeric(enriched_funcs[[balance]][func_id, "P_corr"]) < func_p_cutoff) {
            
            if(as.numeric(enriched_funcs[[balance]][func_id, "OR"]) > 1) {
              func_summaries$balances[[func_id]]$positive_balances <- c(func_summaries$balances[[func_id]]$positive_balances, balance)
            } else if(as.numeric(enriched_funcs[[balance]][func_id, "OR"]) < 1) {
              func_summaries$balances[[func_id]]$negative_balances <- c(func_summaries$balances[[func_id]]$negative_balances, balance)
            } else {
              print(enriched_funcs[[balance]][func_id, ])
              stop("Significant function but OR for above info not different from 1?!")
            }
          } else {
            func_summaries$balances[[func_id]]$nonenriched_sig_balances <- c(func_summaries$balances[[func_id]]$nonenriched_sig_balances, balance)
          }
        } else {
          # Since this balance was NOT significant then categorize it as either enriched or nonenriched.
          if(as.numeric(enriched_funcs[[balance]][func_id, "P_corr"]) < func_p_cutoff) {
            func_summaries$balances[[func_id]]$enriched_nonsig_balances <- c(func_summaries$balances[[func_id]]$enriched_nonsig_balances, balance)
          } else {
            func_summaries$balances[[func_id]]$nonenriched_nonsig_balances <- c(func_summaries$balances[[func_id]]$nonenriched_nonsig_balances, balance)
          }
        }
      } else {
        if(balance %in% sig_balances) {
          func_summaries$balances[[func_id]]$missing_sig_balances <- c(func_summaries$balances[[func_id]]$missing_sig_balances, balance)
        } else {
          func_summaries$balances[[func_id]]$missing_nonsig_balances <- c(func_summaries$balances[[func_id]]$missing_nonsig_balances, balance)
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

two_group_balance_tree_pipeline <- function(abun,
                                            func,
                                            phylogeny,
                                            taxa,
                                            group1_samples,
                                            group2_samples,
                                            ncores=1,
                                            pseudocount=1,
                                            min_num_tips=10,
                                            min_func_freq=10,
                                            min_prop_nonzero=0.1,
                                            name_threshold=0.9,
                                            balance_p_cutoff = 0.05,
                                            balance_correction = "BY",
                                            function_p_cutoff = 0.05,
                                            function_correction = "none",
                                            pathway2func_map = "/home/gavin/projects/functional_reference_frames/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                            pathway_descrip_infile = "/home/gavin/projects/functional_reference_frames/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv",
                                            func_descrip_infile = "/home/gavin/projects/functional_reference_frames/KEGG_mappings/prepped/KO_descrip_22Aug2019.tsv",
                                            verbose=FALSE) {
  
  if(verbose) { message("Prepping input phylogeny.") }
  phylogeny <- prep_tree(phy=phylogeny, tips2keep=rownames(abun))

  if(verbose) { message("Calculating balances.") }
  calculated_balances <- compute_tree_node_balances(abun=abun, phylogeny=phylogeny, ncores=ncores)
  
  # Determine negligible balances based on set cut-offs.
  negligible_balances <- sapply(names(calculated_balances$balances),
                                function(x) {
                                  num_nonzero <- length(which(calculated_balances$balances[[x]] > 0))
                                  total_samples <- length(calculated_balances$balances[[x]])
                                  lhs_feat_num <- length(calculated_balances$features[[x]]$lhs)
                                  rhs_feat_num <- length(calculated_balances$features[[x]]$rhs)
                                  if((num_nonzero < (min_prop_nonzero * total_samples)) || (lhs_feat_num < min_num_tips) || (rhs_feat_num < min_num_tips)) {
                                    return(TRUE)
                                  } else {
                                    return(FALSE) 
                                  }
                                })
  if(verbose) { message("Identified ", length(negligible_balances), " of ", length(calculated_balances$balances), " nodes as negligible and will not be used for analyses.") }

  calculated_balances$balances_filt <- calculated_balances$balances[which(! negligible_balances)]
  
  wilcox_p <- c()
  
  mean_direction <- c()
  
  for(balance in names(calculated_balances$balances_filt)) {
    balance_wilcox <- wilcox.test(calculated_balances$balances_filt[[balance]][group1_samples],
                                  calculated_balances$balances_filt[[balance]][group2_samples],
                                  exact=FALSE)
    wilcox_p <- c(wilcox_p, balance_wilcox$p.value)
    
    group1_mean <- mean(calculated_balances$balances_filt[[balance]][group1_samples])
    group2_mean <- mean(calculated_balances$balances_filt[[balance]][group2_samples])
    
    if(group1_mean > group2_mean) {
      mean_direction <- c(mean_direction, "group1")
    } else if(group1_mean < group2_mean) {
      mean_direction <- c(mean_direction, "group2")
    } else if(group1_mean == group2_mean) {
      mean_direction <- c(mean_direction, "same")
      warning(paste("The calculated balance means are exactly the same for each group at node ", balance, ", which likely indicates a problem.", sep=""))
    }
  }
  
  names(mean_direction) <- names(calculated_balances$balances_filt)
  
  balance_p_corr <- p.adjust(wilcox_p, balance_correction)
  
  if(verbose) { message("Identifying enriched functions at significant nodes.") }
  all_balances_enriched_funcs <- mclapply(names(calculated_balances$balances_filt),
                                          function(x) {
                                            return(node_func_fisher(node = x,
                                                                    in_tree = phylogeny,
                                                                    in_func = func,
                                                                    higher_group=mean_direction[x],
                                                                    pseudocount=1,
                                                                    multiple_test_corr=function_correction))
                                          },
                                          mc.cores = ncores)
  
  names(all_balances_enriched_funcs) <- names(calculated_balances$balances_filt)
  
  sig_balances <- names(calculated_balances$balances_filt)[which(balance_p_corr < balance_p_cutoff)]
  
  if(verbose) { message("Summarizing significant functions across nodes.") }
  func_summaries <- summarize_balance_enrichment(all_balances_enriched_funcs, sig_balances, function_p_cutoff)
  
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
                            "num_sig_balances_nonenrich",
                            "num_sig_balances_pos_enrich",
                            "num_sig_balances_neg_enrich",
                            "num_sig_balances_not_present",
                            "num_nonsig_balances_nonenrich",
                            "num_nonsig_balances_enrich",
                            "num_nonsig_balances_not_present",
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
  
  for(func_id in all_func_id) {
    
    summary_df[func_id, c("num_sig_balances_nonenrich",
                          "num_sig_balances_pos_enrich",
                          "num_sig_balances_neg_enrich",
                          "num_sig_balances_not_present",
                          "num_nonsig_balances_nonenrich",
                          "num_nonsig_balances_enrich",
                          "num_nonsig_balances_not_present")] <- c(length(func_summaries$balances[[func_id]]$nonenriched_sig_balances),
                                                                    length(func_summaries$balances[[func_id]]$positive_balances),
                                                                    length(func_summaries$balances[[func_id]]$negative_balances),
                                                                    length(func_summaries$balances[[func_id]]$missing_sig_balances),
                                                                    length(func_summaries$balances[[func_id]]$nonenriched_nonsig_balances),
                                                                    length(func_summaries$balances[[func_id]]$enriched_nonsig_balances),
                                                                    length(func_summaries$balances[[func_id]]$missing_nonsig_balances))

    phylogeny_node_dists <- dist.nodes(phylogeny)
    
    all_nodes_present <- c(func_summaries$balances[[func_id]]$nonenriched_sig_balances,
                           func_summaries$balances[[func_id]]$positive_balances,
                           func_summaries$balances[[func_id]]$negative_balances,
                           func_summaries$balances[[func_id]]$nonenriched_nonsig_balances,
                           func_summaries$balances[[func_id]]$enriched_nonsig_balances)
    
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
                                                                                          node_labels = func_summaries$balances[[func_id]]$negative_balances),
                                                                  internode_mean_max_dist(phy = phylogeny, dist_matrix = phylogeny_node_dists,
                                                                                          node_labels = func_summaries$balances[[func_id]]$positive_balances))
    
  }
  
  # For each significant balance, test for significantly enriched pathways as well.
  overenriched_balance_funcs <- c()
  underenriched_balance_funcs <- c()
  
  sig_balances_enriched_funcs <- all_balances_enriched_funcs[sig_balances]
  
  for(balance in names(sig_balances_enriched_funcs)) {
    
    balance_sig_func_df <- sig_balances_enriched_funcs[[balance]]
    
    overenriched_balance_funcs <- c(overenriched_balance_funcs, rownames(balance_sig_func_df)[which(balance_sig_func_df$OR > 1 & balance_sig_func_df$P_corr < function_p_cutoff)])
    underenriched_balance_funcs <- c(underenriched_balance_funcs, rownames(balance_sig_func_df)[which(balance_sig_func_df$OR < 1 & balance_sig_func_df$P_corr < function_p_cutoff)])
  }
  
  enriched_balance_funcs <- list(over=overenriched_balance_funcs, under=underenriched_balance_funcs)
  
  if(verbose) { message("Identifying enriched pathways.") }
  sig_balances_enriched_pathways <- identify_enriched_pathways(func=func, func_subsets_to_test=enriched_balance_funcs,
                                                               pathway2func_map=pathway2func_map, pathway_descrip_infile = pathway_descrip_infile)
  
  # Also return list of input parameters.
  input_param <- list(group1_samples=group1_samples,
                      group2_samples=group2_samples,
                      ncores=ncores,
                      pseudocount=pseudocount,
                      min_num_tips=min_num_tips, 
                      min_func_freq=min_func_freq,
                      min_prop_nonzero=min_prop_nonzero,
                      name_threshold=name_threshold,
                      balance_p_cutoff = balance_p_cutoff,
                      balance_correction = balance_correction,
                      function_p_cutoff = function_p_cutoff,
                      function_correction = function_correction,
                      pathway2func_map = pathway2func_map,
                      pathway_descrip_infile = pathway_descrip_infile,
                      func_descrip_infile = func_descrip_infile)
  
  return(list(balances_info=calculated_balances,
              sig_balances=sig_balances,
              pathways=sig_balances_enriched_pathways,
              funcs_per_balance=sig_balances_enriched_funcs,
              funcs=enriched_balance_funcs,
              df=summary_df,
              out_list=func_summaries,
              tree=phylogeny,
              input_param=input_param))
  
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

