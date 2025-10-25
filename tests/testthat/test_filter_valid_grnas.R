test_that("filter_valid_grnas() returns GRanges subset", {
  skip_on_cran()
  expect_error(filter_valid_grnas(mock_gr, genome=NULL, species="test"))
})