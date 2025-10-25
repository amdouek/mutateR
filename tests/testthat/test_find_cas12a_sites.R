test_that("find_cas12a_sites() returns GRanges with 34bp context", {
  skip_on_cran()
  library(BSgenome.Hsapiens.UCSC.hg38)
  exon_gr <- get_exon_structures("ENST00000269305", "hsapiens", "GRanges")
  hits <- find_cas12a_sites(exon_gr[1], BSgenome.Hsapiens.UCSC.hg38)
  expect_s4_class(hits, "GRanges")
  expect_equal(unique(nchar(mcols(hits)$sequence_context)), 34)
})