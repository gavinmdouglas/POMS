
ex_tree_w_label <- ex_tree
ex_tree_wo_label <- ex_tree
ex_tree_wo_label$node.label <- NULL

expected_df <- data.frame(num_FSNs = c(4, 5),
                          num_FSNs_group1_enrich = c(2, 5),
                          num_FSNs_group2_enrich = c(1, 0),
                          num_FSNs_at_nonBSNs = c(1, 0),
                          multinomial_p = c(1, 0.0272),
                          multinomial_corr = c(1, 0.0544))
rownames(expected_df) <- c("K07106", "K02036")

ex_balances <- compute_node_balances(tree = ex_tree_w_label,
                                     abun_table = ex_taxa_abun,
                                     min_num_tips = 4,
                                     ncores = 1,
                                     pseudocount = 1)


test_that("two-group pipeline produces expected basic output with example files", {
  
  ex_basic_output <- POMS_pipeline(abun = ex_taxa_abun,
                                    func = ex_func,
                                    tree = ex_tree_wo_label,
                                    group1_samples = ex_group1,
                                    group2_samples = ex_group2,
                                    ncores = 1,
                                    min_num_tips = 4,
                                    multinomial_min_FSNs = 3,
                                    min_func_instances = 0,
                                   derep_nodes = FALSE)

  expect_equal(ex_basic_output$results, expected_df)
})


test_that("two-group pipeline produces expected basic output with example files and manually specified BSNs and balances.", {
  
  ex_basic_output <- POMS_pipeline(abun = ex_taxa_abun,
                                   func = ex_func,
                                   tree = ex_tree_w_label,
                                   group1_samples = ex_group1,
                                   group2_samples = ex_group2,
                                   manual_BSNs = c("n1", "n2", "n13", "n17", "n18", "n20", "n32", "n46"),
                                   manual_balances = ex_balances$balances,
                                   ncores = 1,
                                   min_num_tips = 4,
                                   multinomial_min_FSNs = 3,
                                   min_func_instances = 0,
                                   derep_nodes = FALSE)
  
  expect_equal(ex_basic_output$results, expected_df)
})


test_that("two-group pipeline produces expected basic output with example files and manually specified BSNs, balances, and BSN directions.", {
  
  manual_BSN_dir_set <- c("group1", "group2", "group1", "group1", "group2", "group1", "group2", "group1")
  names(manual_BSN_dir_set) <- c("n1", "n2", "n13", "n17", "n18", "n20", "n32", "n46")
  
  ex_basic_output <- POMS_pipeline(abun = ex_taxa_abun,
                                   func = ex_func,
                                   tree = ex_tree_w_label,
                                   manual_BSNs = c("n1", "n2", "n13", "n17", "n18", "n20", "n32", "n46"),
                                   manual_balances = ex_balances$balances,
                                   manual_BSN_dir = manual_BSN_dir_set,
                                   ncores = 1,
                                   min_num_tips = 4,
                                   multinomial_min_FSNs = 3,
                                   min_func_instances = 0,
                                   derep_nodes = FALSE)
  
  expect_equal(ex_basic_output$results, expected_df)
})


test_that("correct error occurs when significant nodes are not subset of tested nodes.", {
  
  expect_error(object = POMS_pipeline(abun = ex_taxa_abun,
                                      func = ex_func,
                                      tree = ex_tree_w_label,
                                      group1_samples = ex_group1,
                                      group2_samples = ex_group2,
                                      ncores = 1,
                                      min_num_tips = 4,
                                      multinomial_min_FSNs = 3,
                                      min_func_instances = 0,
                                      manual_BSNs = c("test1", "test2"),
                                      manual_balances = list("a"=as.numeric(), "b"=as.numeric(), "c"=as.numeric()),
                                      derep_nodes = FALSE),
               regexp = "not all nodes in manual_BSNs vector are present in manual_balances object")
})



test_that("correct error occurs when some node labels in balances input are not found in tree.", {
  
  ex_balances_prepped <- compute_node_balances(tree = ex_tree_w_label,
                                               abun_table = ex_taxa_abun,
                                               min_num_tips=5,
                                               ncores=1,
                                               pseudocount=1)
  
  names(ex_balances_prepped$balances)[1] <- "test"
  
  expect_error(object = POMS_pipeline(abun = ex_taxa_abun,
                                      func = ex_func,
                                      tree = ex_tree_w_label,
                                      group1_samples = ex_group1,
                                      group2_samples = ex_group2,
                                      ncores = 1,
                                      min_num_tips = 4,
                                      multinomial_min_FSNs = 3,
                                      min_func_instances = 0,
                                      manual_BSNs = c("test", "n2"),
                                      manual_balances = ex_balances_prepped$balances,
                                      derep_nodes = FALSE),
                                      regexp = "Stopping - some balance labels \\(in manually input manual_balances argument\\) are missing from tree node labels.")
})



test_that("correct error occurs when input tree is missing node labels when manual BSNs are specified.", {
  
  expect_error(object = POMS_pipeline(abun = ex_taxa_abun,
                                      func = ex_func,
                                      tree = ex_tree_wo_label,
                                      group1_samples = ex_group1,
                                      group2_samples = ex_group2,
                                      ncores = 1,
                                      min_num_tips = 4,
                                      multinomial_min_FSNs = 3,
                                      min_func_instances = 0,
                                      manual_BSNs = c("a", "b"),
                                      manual_balances = list("a"=as.numeric(), "b"=as.numeric(), "c"=as.numeric()),
                                      derep_nodes = FALSE),
                                      regexp = "Stopping - node labels must be present in tree if manual balances \\(i.e., manual_balances argument\\) are specificed.")
})


test_that("prep_func_node_info works as expected with pipeline output", {
  
  ex_basic_output <- POMS_pipeline(abun = ex_taxa_abun,
                                   func = ex_func,
                                   tree = ex_tree_wo_label,
                                   group1_samples = ex_group1,
                                   group2_samples = ex_group2,
                                   ncores = 1,
                                   min_num_tips = 4,
                                   multinomial_min_FSNs = 3,
                                   min_func_instances = 0,
                                   derep_nodes = FALSE)
  
  K07106_node_info <- prep_func_node_info(POMS_output = ex_basic_output,
                                          func_id = "K07106",
                                          taxa_table = NULL)
  
  expected_node_info_wo_tree <- list(FSNs_group1_enrich_i = c(78, 80),
                                     FSNs_group2_enrich_i = 61,
                                     FSNs_at_nonBSNs_i = 93,
                                     all_FSNs_i = c(78, 80, 61, 93),
                                     all_BSNs_i = c(61, 62, 73, 77, 78, 80, 92, 106),
                                     BSNs_at_nonFSNs_i = c(62, 73, 77, 92, 106),
                                     tested_nodes_i = c(61, 62, 73, 77, 78, 80, 85, 92, 93, 106))
  
  expect_equal(K07106_node_info[names(expected_node_info_wo_tree)], expected_node_info_wo_tree)
})

