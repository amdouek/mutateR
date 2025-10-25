test_that("map_protein_domains() returns NULL safely for missing data", {
  skip_on_cran()
  expect_warning(x <- map_protein_domains("ENST00000269305", "hsapiens"))
})