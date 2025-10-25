test_that("check_frameshift_ptc() detects frameshift correctly", {
  df <- data.frame(rank=1:5, exon_cds_length=c(150,150,90,120,120))
  res <- check_frameshift_ptc(df, 1, 4)
  expect_named(res, c("frameshift","ptc_flag","terminal_exon_case"))
  expect_type(res$frameshift, "logical")
})