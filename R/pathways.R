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


parse_sig_pathway_names <- function(pathway_list) {
  
  out_names <- c()
  
  for(pathway in names(pathway_list)) {
    out_names <- c(out_names, paste(pathway, pathway_list[[pathway]]$descrip, sep=" - "))
  }
  return(out_names)
}


gene_set_enriched_pathways <- function(up_genes,
                                       down_genes,
                                       gene_background,
                                       p_corr_method="BY",
                                       corr_P_cutoff=0.05,
                                       pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                       pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv") {
  
  
  focal_genes <- c(up_genes, down_genes)
  
  if(! file.exists(pathway2func_map)) { stop("Stopping - file corresponding to pathway2func_map argument not found.") } 
  if(! file.exists(pathway_descrip_infile)) { stop("Stopping - file corresponding to pathway_descrip_infile argument not found.") } 
  if(length(which(! focal_genes %in% gene_background)) > 0) { stop("Stopping - all genes in focal set must also be in gene background.") } 
  
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
    
    if(length(which(! pathway_funcs %in% gene_background)) > 0) {
      pathway_funcs <- pathway_funcs[-which(! pathway_funcs %in% gene_background)]
    }
    
    pathway2funcsets[[pathway]] <- pathway_funcs
  }
  
  # Get background of how many functions in this dataset are in each pathway.
  total_funcs <- length(gene_background)
  background_pathway_func_counts <- c()
  
  for(pathway in names(pathway2funcsets)) {
    func_set <- pathway2funcsets[[pathway]]
    func_set <- func_set[which(func_set %in% gene_background)]
    
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
  
  func_subsets_to_test <- list(up=up_genes,
                               down=down_genes)
  
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
