test_that("find_cas9_sites() returns GRanges with required fields", {
  skip_on_cran()
  library(BSgenome.Hsapiens.UCSC.hg38)
  exon_gr <- get_exon_structures("ENST00000269305", "hsapiens", "GRanges")
  hits <- find_cas9_sites(exon_gr[1], BSgenome.Hsapiens.UCSC.hg38)
  expect_s4_class(hits, "GRanges")
  expect_true(all(c("protospacer_sequence","pam_sequence") %in% names(mcols(hits))))
})

test_that("find_cas9_sites() reports correct genomic coords on both strands", {
  # Regression test for the minus-strand coordinate bug fixed 2026-05-11:
  # Previously the minus-strand GRanges range was placed ~22 bp downstream
  # of the actual protospacer and cut_site was ~8 bp off. The invariant
  # checked here — that getSeq(genome, hit) returns the recorded
  # target_sequence — is strand-aware (BSgenome respects the GRanges strand)
  # and would have caught the bug.
  skip_on_cran()
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38)
  exon_gr <- get_exon_structures("ENST00000269305", "hsapiens", "GRanges")
  hits <- find_cas9_sites(exon_gr[1], BSgenome.Hsapiens.UCSC.hg38)
  expect_true(length(hits) > 0)

  # Each hit's GRanges + strand should fetch the same target_sequence (proto+PAM)
  fetched <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, hits))
  recorded <- as.character(mcols(hits)$target_sequence)
  expect_equal(fetched, recorded)

  # cut_site must fall inside the reported genomic range for every hit
  expect_true(all(BiocGenerics::start(hits) <= mcols(hits)$cut_site &
                    mcols(hits)$cut_site <= BiocGenerics::end(hits)))

  # Both strands should be represented in a non-trivial exon
  expect_true(any(as.character(strand(hits)) == "+"))
  expect_true(any(as.character(strand(hits)) == "-"))
})