
library(POMS)


ex_taxa_abun <- read.table("../../example_files/ex_taxa_abun.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
ex_taxa_abun <- ex_taxa_abun[, c("ERR321132", "SRR5127690", "SRR2992922", "ERR321066")]

ex_tree <- ape::read.tree("../../example_files/ex_tree.newick")

test_features1 <- c("SRR769509_bin.3", "21673_4_5", "SRR3992981_bin.16")
test_features2 <- c("ERR1620320_bin.20", "ERR1293861_bin.21", "SRR4408017_bin.22", "SRR2155382_bin.36", "ERR1620320_bin.43")

test_features_w_missing <- c("ERR1620320_bin.20", "ERR1293861_bin.21", "not_present")
test_features_w_overlapping <- c("ERR1620320_bin.20", "ERR1293861_bin.21", "21673_4_5")

ex_taxa_labels <- read.table("../../example_files/ex_taxa_labels.tsv.gz",
                             header = TRUE,
                             sep = "\t",
                             stringsAsFactors = FALSE,
                             row.names = 1)

test_that("Expected error returned when features not present in table are input.", {
  expect_error(object = abun_isometric_log_ratios(abun_table = ex_taxa_abun,
                                                  set1_features = test_features1,
                                                  set2_features = test_features_w_missing),
               regexp = "Stopping - at least one feature in the specified sets is not present as a row name in the abundance table.")
})


test_that("Expected error returned when features intersecting between sets.", {
  expect_error(object = abun_isometric_log_ratios(abun_table = ex_taxa_abun,
                                                  set1_features = test_features1,
                                                  set2_features = test_features_w_overlapping),
               regexp = "Stopping - at least one feature overlaps between the input sets.")
})


test_that("Check that error related to presence of 0's occurs in absence of pseudocount", {
  expect_error(object = abun_isometric_log_ratios(abun_table = ex_taxa_abun,
                                                  set1_features = test_features1,
                                                  set2_features = test_features2,
                                                  pseudocount = 0),
               regexp = "At least one 0 is present in the abundance table, which means that at least some isometric log ratios cannot be computed.")
})


test_that("Check that ILR values make sense for test case.", {
  
  expected_output <- c(1.622685, 1.279992, 0.000000, -0.012730)
  names(expected_output) <- c("ERR321132", "SRR5127690", "SRR2992922", "ERR321066")
  
  expect_equal(object = abun_isometric_log_ratios(abun_table = ex_taxa_abun,
                                                  set1_features = test_features1,
                                                  set2_features = test_features2,
                                                  pseudocount = 1),
               expected = expected_output,
               tolerance = 0.000001)
})



test_that("compute_node_balances error when no node labels.", {
  
  ex_tree$node.label <- NULL
  
  expect_error(object = compute_node_balances(phylogeny = ex_tree,
                                              abun_table = ex_taxa_abun,
                                              min_num_tips=5,
                                              ncores=1,
                                              pseudocount=1),
               regexp = "Stopping - input tree does not have any node labels.")

})


test_that("compute_node_balances error when tip not found as row name.", {
  
  ex_tree$tip.label[1] <- "test"
  
  expect_error(object = compute_node_balances(phylogeny = ex_tree,
                                              abun_table = ex_taxa_abun,
                                              min_num_tips=5,
                                              ncores=1,
                                              pseudocount=1),
               regexp = "Stopping - not all tips are found as row names in the abundance table.")
  
})


test_that("compute_node_balances error when node name in subset_to_test not found in tree.", {
  
  expect_error(object = compute_node_balances(phylogeny = ex_tree,
                                              abun_table = ex_taxa_abun,
                                              min_num_tips=5,
                                              ncores=1,
                                              pseudocount=1,
                                              subset_to_test = c("n1", "n2", "n1000")),
               regexp = "Stopping - some labels in subset_to_test do not match node labels in the tree.")
  
})


test_that("compute_node_balances check balances at one node with all default settings.", {
  
  balances_out <- compute_node_balances(phylogeny = ex_tree,
                                        abun_table = ex_taxa_abun,
                                        min_num_tips=5,
                                        ncores=1,
                                        pseudocount=1)
  
  expect_equal(as.numeric(balances_out$balances$n2),
               c(0.4557199, 0.1411240, 0.3392601, 0.1679282),
               tolerance = 0.0000001)
  
})


test_that("compute_node_balances change min number of tips to make sure that's working.", {
  
  balances_out <- compute_node_balances(phylogeny = ex_tree,
                                        abun_table = ex_taxa_abun,
                                        min_num_tips=2,
                                        ncores=1,
                                        pseudocount=1)
  
  expect_equal(as.numeric(balances_out$balances$n47),
               c(0.2362104, -0.5592617, -0.7025587, 0.0000000),
               tolerance = 0.0000001)
  
})


test_that("compute_node_balances try subset_to_test with nodes present and no nodes that pass the min_num_tips", {
  
  expect_error(object = compute_node_balances(phylogeny = ex_tree,
                                              abun_table = ex_taxa_abun,
                                              min_num_tips=5,
                                              ncores=1,
                                              pseudocount=1,
                                              subset_to_test = c("n3", "n8", "n13")),
               regexp = "Stopping - no non-negligible nodes remain after filtering based on mininum number of tips of left and right-hand side of each node.")
  
})


test_that("compute_node_balances try altering pseudocount to make sure that parameter is working as expected.", {
  
  balances_out <- compute_node_balances(phylogeny = ex_tree,
                                        abun_table = ex_taxa_abun,
                                        min_num_tips=5,
                                        ncores=1,
                                        pseudocount=0.1)
  
  expect_equal(as.numeric(balances_out$balances$n2),
               c(0.6441389, 0.2443771, 0.9872499, 0.3890996),
               tolerance = 0.0000001)
  
})


test_that("Check taxa on each side of the node - combined labels and threshold of 0.51", {
 
  exp_output <- c("Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Porphyromonadaceae (Family)",
                  "Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Porphyromonadaceae; NA (Genus)")

  obs_output <- node_taxa(in_tree = ex_tree, taxon_labels = ex_taxa_labels,
                          node_label = "n1", combine_labels = TRUE, threshold = 0.51)
  
  expect_equal(obs_output, exp_output)

})


test_that("Check taxa on each side of the node - simple (non-combined) labels and threshold of 0.95", {
  
  exp_output <- c("Porphyromonadaceae (Family)", "Bacteroidales (Order)")
  
  obs_output <- node_taxa(in_tree = ex_tree, taxon_labels = ex_taxa_labels,
                          node_label = "n1", combine_labels = FALSE, threshold = 0.95)
  
  expect_equal(obs_output, exp_output)
  
})
