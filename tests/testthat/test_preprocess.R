# Tests for small preprocessing functions.

library(POMS)

test_input <- data.frame(matrix(c(0, 1, 30, 0.1, 10, 0, 0, 0, 0, 1, 0, 10), nrow = 4, ncol = 3))
rownames(test_input) <- c("row1", "row2", "row3", "row4")
colnames(test_input) <- c("colA", "colB", "colC")


test_that("filter_rare_table_cols is working as expected for one example (based on min *count* and not dropping missing rows)", {
  
  expect_equal(filter_rare_table_cols(in_tab=test_input, min_nonzero_count=2, min_nonzero_prop=0.01, drop_missing_rows = FALSE, verbose=FALSE),
               test_input[, -2])
})
 

test_that("filter_rare_table_cols is working as expected for one example (based on min *proportion* and not dropping missing rows)", {
   
   expect_equal(filter_rare_table_cols(in_tab=test_input, min_nonzero_count=1, min_nonzero_prop=0.26, drop_missing_rows=FALSE, verbose=FALSE),
                test_input[, -2])
})


test_that("filter_rare_table_cols is working as expected for one example (based on min *count* and dropping missing rows)", {
  
  expect_equal(filter_rare_table_cols(in_tab=test_input, min_nonzero_count=2, min_nonzero_prop=0.01, drop_missing_rows = TRUE, verbose=FALSE),
               test_input[-1, -2])
})


test_that("subset_by_col_and_filt is working as expected for one example", {
  
  expect_equal(subset_by_col_and_filt(in_tab=test_input, col2keep=c("colA", "colC"), verbose=FALSE),
               test_input[-1, -2])
})
