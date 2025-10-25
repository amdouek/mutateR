test_that("score_grnas() adds gc and ontarget_score columns", {
  seqs <- Biostrings::DNAStringSet(rep("AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", 3))
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:3, width=30))
  mcols(gr)$sequence_context <- as.character(seqs)
  res <- score_grnas(gr, method="ruleset1")
  expect_true(all(c("gc","ontarget_score") %in% names(mcols(res))))
})