test_that("find_cas9_sites() returns GRanges with required fields", {
  skip_on_cran()
  library(BSgenome.Hsapiens.UCSC.hg38)
  exon_gr <- get_exon_structures("ENST00000269305", "hsapiens", "GRanges")
  hits <- find_cas9_sites(exon_gr[1], BSgenome.Hsapiens.UCSC.hg38)
  expect_s4_class(hits, "GRanges")
  expect_true(all(c("protospacer_sequence","pam_sequence") %in% names(mcols(hits))))
})