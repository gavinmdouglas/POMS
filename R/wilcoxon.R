#' @export
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

#' @export
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

