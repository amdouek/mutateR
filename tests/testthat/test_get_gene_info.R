test_that("get_gene_info() returns canonical and all data frames", {
  skip_on_cran()
  expect_error(get_gene_info(), "species")
  result <- get_gene_info("TP53", "hsapiens")
  expect_type(result, "list")
  expect_true(all(c("canonical","all") %in% names(result)))
})