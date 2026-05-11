test_that("nuclease_spec() builds a valid SauCas9 spec", {
  s <- nuclease_spec(
    name               = "SauCas9",
    pam                = "NNGRRT",
    pam_side           = "3prime",
    protospacer_length = 21L,
    cut_offset_top     = -3L,
    cut_offset_bottom  = -3L,
    activity           = "cut",
    grna_architecture  = "crRNA + tracr",
    source             = "CasPEDIA:1.1.1"
  )
  expect_s3_class(s, "NucleaseSpec")
  expect_equal(s$pam, "NNGRRT")
  expect_equal(s$pam_length, 6L)
  expect_equal(s$pam_side, "3prime")
  expect_equal(s$protospacer_length, 21L)
  expect_false(s$is_canonical)
})

test_that("nuclease_spec() builds a valid staggered Cas12b spec", {
  s <- nuclease_spec(
    name               = "AacCas12b",
    pam                = "TTN",
    pam_side           = "5prime",
    protospacer_length = 20L,
    cut_offset_top     = -5L,
    cut_offset_bottom  = -9L
  )
  expect_equal(s$cut_offset_top, -5L)
  expect_equal(s$cut_offset_bottom, -9L)
  expect_false(s$is_canonical)
})

test_that("nuclease_spec() builds a valid nickase spec (Cas12i)", {
  s <- nuclease_spec(
    name               = "Cas12i1",
    pam                = "TTN",
    pam_side           = "5prime",
    protospacer_length = 28L,
    cut_offset_top     = -5L,
    cut_offset_bottom  = -5L,
    activity           = "nick"
  )
  expect_equal(s$activity, "nick")
})

test_that("nuclease_spec() rejects RNA-targeting effectors", {
  expect_error(
    nuclease_spec(
      name               = "LbuCas13a",
      pam                = "N",
      pam_side           = "5prime",
      protospacer_length = 20L,
      cut_offset_top     = 0L,
      cut_offset_bottom  = 0L,
      target_type        = "RNA",
      activity           = "cut"
    ),
    regexp = "dsDNA-targeting"
  )
})

test_that("nuclease_spec() rejects binding-only effectors", {
  expect_error(
    nuclease_spec(
      name               = "Cas12k",
      pam                = "GTN",
      pam_side           = "5prime",
      protospacer_length = 20L,
      cut_offset_top     = 0L,
      cut_offset_bottom  = 0L,
      activity           = "bind"
    ),
    regexp = "'cut' or 'nick'"
  )
})

test_that("nuclease_spec() rejects unknown-activity effectors", {
  expect_error(
    nuclease_spec(
      name               = "Cas12U2",
      pam                = "N",
      pam_side           = "5prime",
      protospacer_length = 20L,
      cut_offset_top     = 0L,
      cut_offset_bottom  = 0L,
      activity           = "unknown"
    ),
    regexp = "'cut' or 'nick'"
  )
})

test_that("nuclease_spec() rejects malformed PAM strings", {
  expect_error(
    nuclease_spec(
      name               = "BadPAM",
      pam                = "XYZ",
      pam_side           = "5prime",
      protospacer_length = 20L,
      cut_offset_top     = 0L,
      cut_offset_bottom  = 0L
    ),
    regexp = "IUPAC"
  )
})

test_that("nuclease_spec() rejects out-of-range protospacer length", {
  expect_error(
    nuclease_spec(
      name               = "TooShort",
      pam                = "NGG",
      pam_side           = "3prime",
      protospacer_length = 5L,
      cut_offset_top     = 0L,
      cut_offset_bottom  = 0L
    ),
    regexp = "protospacer_length"
  )
})

test_that("nuclease_spec() rejects pam_is_pfs = TRUE (reserved)", {
  expect_error(
    nuclease_spec(
      name               = "Cas13b-like",
      pam                = "NAN",
      pam_side           = "3prime",
      protospacer_length = 30L,
      cut_offset_top     = 0L,
      cut_offset_bottom  = 0L,
      pam_is_pfs         = TRUE
    ),
    regexp = "pam_is_pfs"
  )
})

test_that("resolve_nuclease() returns canonical specs for legacy strings", {
  for (key in c("Cas9", "Cas12a", "enCas12a")) {
    s <- mutateR:::resolve_nuclease(key)
    expect_s3_class(s, "NucleaseSpec")
    expect_true(s$is_canonical)
    expect_equal(s$canonical_key, key)
  }
})

test_that("resolve_nuclease() passes through NucleaseSpec inputs", {
  custom <- nuclease_spec(
    name               = "TestSpec",
    pam                = "NGG",
    pam_side           = "3prime",
    protospacer_length = 20L,
    cut_offset_top     = -3L,
    cut_offset_bottom  = -3L
  )
  expect_identical(mutateR:::resolve_nuclease(custom), custom)
})

test_that("resolve_nuclease() rejects unknown strings", {
  expect_error(mutateR:::resolve_nuclease("BoringNuclease"))
})

test_that("print.NucleaseSpec runs without error", {
  s <- nuclease_spec(
    name               = "AacCas12b",
    pam                = "TTN",
    pam_side           = "5prime",
    protospacer_length = 20L,
    cut_offset_top     = -5L,
    cut_offset_bottom  = -9L
  )
  expect_output(print(s), "NucleaseSpec")
  expect_output(print(s), "staggered")
})
