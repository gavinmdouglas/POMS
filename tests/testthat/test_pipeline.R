# Test that the overall pipeline function is producing the expected output.

library(POMS)

ex1_taxa_abun <- read.table("../../example_files/ex_taxa_abun.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
ex1_func <- read.table("../../example_files/ex_func.tsv.gz", header = TRUE, sep = "\t", row.names = 1)

ex1_group1 <- read.table("../../example_files/ex_group1.txt.gz", stringsAsFactors = FALSE)$V1
ex1_group2 <- read.table("../../example_files/ex_group2.txt.gz", stringsAsFactors = FALSE)$V1

ex1_tree <- ape::read.tree("../../example_files/ex_tree.newick")
ex1_tree_w_label <- ex1_tree
ex1_tree$node.label <- NULL

# Example of how to run main POMS function. 
test_that("Two-group pipeline produces expected basic output with ex1 files", {
 
  expected_df <- data.frame(num_FSNs = c(4, 5),
                            num_FSNs_group1_enrich = c(2, 5),
                            num_FSNs_group2_enrich = c(1, 0),
                            num_FSNs_at_nonBSNs = c(1, 0),
                            multinomial_p = c(1, 0.0272),
                            multinomial_corr = c(1, 0.0544))
  rownames(expected_df) <- c("K07106", "K02036")
  
  ex1_basic_output <- POMS_pipeline(abun = ex1_taxa_abun,
                                    func = ex1_func,
                                    phylogeny = ex1_tree,
                                    group1_samples = ex1_group1,
                                    group2_samples = ex1_group2,
                                    ncores = 1,
                                    min_num_tips = 4,
                                    multinomial_min_FSNs = 3,
                                    min_func_instances = 0)

  expect_equal(ex1_basic_output$summary, expected_df)
})
