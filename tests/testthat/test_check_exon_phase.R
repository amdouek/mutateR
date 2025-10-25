test_that("check_exon_phase() identifies phase-compatible exon pairs", {
  df <- data.frame(rank = 1:3,
                   start_phase = c(0,2,0),
                   end_phase =   c(2,0,0))
  res <- check_exon_phase(df)
  expect_true(is.data.frame(res))
  expect_true(all(c("exon_5p","exon_3p","compatible") %in% names(res)))
})