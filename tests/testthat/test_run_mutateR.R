test_that("run_mutateR() completes pipeline without errors (Cas9)", {
  skip_on_cran()
  library(BSgenome.Hsapiens.UCSC.hg38)
  res <- run_mutateR("TP53", "hsapiens",
                     genome = BSgenome.Hsapiens.UCSC.hg38,
                     nuclease = "Cas9",
                     top_n = 1)
  expect_type(res, "list")
  expect_true(all(c("pairs","plot") %in% names(res)))
})