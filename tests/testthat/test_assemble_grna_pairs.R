test_that("assemble_grna_pairs() returns a data.frame", {
  skip_on_cran()
  expect_error(assemble_grna_pairs(mock_gr, mock_gr, "txid", "species"))
})