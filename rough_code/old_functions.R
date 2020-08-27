library(grid)
library(VennDiagram)
library(corncob)
library(ALDEx2)

identify_enriched_pathways <- function(func, func_subsets_to_test, pathway2func_map, pathway_descrip_infile, p_corr_method="BH") {
  
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
  
  # Get background of how many functions across genomes in this dataset are in each pathway.
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
    
    pathway_fisher_tests_corr <- p.adjust(pathway_fisher_tests, p_corr_method)
    
    for(i in which(pathway_fisher_tests_corr < 0.05)) {
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
