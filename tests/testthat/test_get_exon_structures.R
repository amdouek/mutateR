test_that("get_exon_structures() properly handles transcript input", {
  skip_on_cran()
  tx_info <- get_gene_info("TP53", "hsapiens")
  tx_id <- tx_info$canonical$ensembl_transcript_id[1]
  exons <- get_exon_structures(tx_id, "hsapiens")
  expect_s3_class(exons, "data.frame")
  expect_true(all(c("start_phase","end_phase") %in% names(exons)))
})