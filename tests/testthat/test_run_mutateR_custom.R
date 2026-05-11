test_that("run_mutateR_custom() returns expected structure (smoke test)", {
  skip_on_cran()
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38)

  sau <- nuclease_spec(
    name               = "SauCas9",
    pam                = "NNGRRT",
    pam_side           = "3prime",
    protospacer_length = 21L,
    cut_offset_top     = -3L,
    cut_offset_bottom  = -3L,
    source             = "CasPEDIA:1.1.1"
  )

  result <- run_mutateR_custom(
    gene_id        = "TP53",
    species        = "hsapiens",
    genome         = BSgenome.Hsapiens.UCSC.hg38,
    nuclease_spec  = sau,
    assemble_pairs = TRUE,
    design_primers = FALSE,
    plot           = FALSE,
    quiet          = TRUE
  )

  expect_type(result, "list")
  expect_true(all(c("gene_id", "transcript_id", "nuclease",
                    "nuclease_spec", "exons", "sites", "pairs") %in%
                  names(result)))
  expect_s3_class(result$nuclease_spec, "NucleaseSpec")
  expect_equal(result$nuclease, "SauCas9")
  if (!is.null(result$sites)) {
    expect_s4_class(result$sites, "GRanges")
    # No scoring should have been performed
    expect_true("ontarget_score" %in%
                  names(GenomicRanges::mcols(result$sites)))
    expect_true(all(is.na(GenomicRanges::mcols(result$sites)$ontarget_score)))
    expect_equal(unique(GenomicRanges::mcols(result$sites)$scoring_method),
                 "none")
  }
})

test_that("run_mutateR() accepts a custom NucleaseSpec and skips scoring", {
  skip_on_cran()
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38)

  cas12i <- nuclease_spec(
    name               = "Cas12i1",
    pam                = "TTN",
    pam_side           = "5prime",
    protospacer_length = 28L,
    cut_offset_top     = -5L,
    cut_offset_bottom  = -5L,
    activity           = "nick",
    source             = "CasPEDIA:4.4.3"
  )

  result <- run_mutateR(
    gene_id   = "TP53",
    species   = "hsapiens",
    genome    = BSgenome.Hsapiens.UCSC.hg38,
    nuclease  = cas12i,
    offtarget = FALSE,
    design_primers = FALSE,
    quiet     = TRUE
  )
  expect_equal(result$nuclease, "Cas12i1")
  expect_equal(result$score_method, "none")
  expect_s3_class(result$nuclease_spec, "NucleaseSpec")
  if (!is.null(result$scored_grnas)) {
    expect_true(all(is.na(GenomicRanges::mcols(result$scored_grnas)$ontarget_score)))
  }
})
