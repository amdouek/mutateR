#' @title Bowtie off-target search backend
#'
#' @description Internal functions for off-target searching via crisprBowtie.
#' Provides fast mismatch-only genome-wide alignment (typically <1 min for
#' ~400 guides at 3 mismatches against the human genome) as the default
#' off-target backend. Results are normalised to the same hit schema used
#' by the CRISPRitz backend for downstream scoring by
#' \code{\link{score_and_aggregate}}.
#'
#' @keywords internal
#' @name ot_bowtie
NULL


#' Map mutateR nuclease names to crisprBase CrisprNuclease objects
#'
#' Returns the appropriate pre-built \code{CrisprNuclease} object from
#' \pkg{crisprBase}. Falls back to manual construction for \code{enAsCas12a}
#' if not available in older crisprBase versions.
#'
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a".
#'
#' @return A \code{CrisprNuclease} object.
#' @keywords internal
get_crisprbase_nuclease <- function(nuclease) {

  if (!requireNamespace("crisprBase", quietly = TRUE)) {
    stop("Package 'crisprBase' is required for the Bowtie off-target backend.\n",
         "Install with: BiocManager::install('crisprBase')")
  }

  nuc_obj <- switch(nuclease,

                    "Cas9" = tryCatch(
                      crisprBase::SpCas9,
                      error = function(e) {
                        stop("Could not load crisprBase::SpCas9. ",
                             "Ensure crisprBase is installed and up to date.\n",
                             "Install with: BiocManager::install('crisprBase')")
                      }
                    ),

                    "Cas12a" = tryCatch(
                      crisprBase::AsCas12a,
                      error = function(e) {
                        stop("Could not load crisprBase::AsCas12a. ",
                             "Ensure crisprBase is installed and up to date.\n",
                             "Install with: BiocManager::install('crisprBase')")
                      }
                    ),

                    "enCas12a" = tryCatch(
                      crisprBase::enAsCas12a,
                      error = function(e) {
                        # Fallback: construct enAsCas12a manually for older crisprBase versions
                        # enAsCas12a uses TTTN PAM (5'), 23-nt spacer, same cut pattern as AsCas12a
                        message("crisprBase::enAsCas12a not found. ",
                                "Constructing enAsCas12a nuclease object manually.")
                        tryCatch({
                          # Start from AsCas12a as template for cut pattern and spacer length
                          base_nuc <- crisprBase::AsCas12a

                          crisprBase::CrisprNuclease(
                            nucleaseName  = "enAsCas12a",
                            targetType    = "DNA",
                            pams          = c("TTTA", "TTTC", "TTTG", "TTTT"),
                            weights       = rep(1, 4),
                            metadata      = list(description = "enAsCas12a (TTTN PAM)"),
                            pam_side      = "5prime",
                            spacer_gap    = crisprBase::spacerGap(base_nuc),
                            spacer_length = crisprBase::spacerLength(base_nuc)
                          )
                        }, error = function(e2) {
                          stop("Failed to construct enAsCas12a nuclease object: ", e2$message, "\n",
                               "Please update crisprBase: BiocManager::install('crisprBase')")
                        })
                      }
                    ),

                    stop("Unsupported nuclease for Bowtie backend: ", nuclease)
  )

  nuc_obj
}


#' Execute Bowtie genome-wide off-target search via crisprBowtie
#'
#' Thin wrapper around \code{crisprBowtie::runCrisprBowtie()} that maps
#' mutateR nuclease types to crisprBase objects and provides sensible
#' defaults for off-target discovery (spacer mode, non-canonical PAMs
#' included). Returns the raw crisprBowtie output data.frame.
#'
#' @param unique_spacers Character vector. Unique protospacer sequences to
#'        search (all must be the same length).
#' @param genome BSgenome object. Required for spacer-mode PAM validation.
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a".
#' @param bowtie_index Character. Path to Bowtie index prefix, as returned
#'        by \code{\link{ensure_bowtie_index}}.
#' @param max_mismatches Integer. Maximum mismatches (0–3). Values > 3 are
#'        clamped with a warning (Bowtie hard limit).
#' @param all_alignments Logical. Return all possible alignments? (default TRUE).
#'        Set to FALSE with \code{n_max_alignments} to cap output volume.
#' @param n_max_alignments Integer. Maximum alignments per spacer when
#'        \code{all_alignments = FALSE} (default 1000).
#' @param canonical_only Logical. Return only hits adjacent to canonical PAM
#'        sequences? (default TRUE). The MIT specificity formula (Hsu et al., 2013)
#'        was calibrated against NGG PAMs only. Set to FALSE for comprehensive
#'        non-canonical PAM analysis; however, this drastically increases
#'        hit volume and reduces specificity scores.
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return A data.frame with columns as returned by
#'   \code{crisprBowtie::runCrisprBowtie()}:
#'   \describe{
#'     \item{spacer}{Input spacer sequence.}
#'     \item{protospacer}{Genomic off-target protospacer sequence.}
#'     \item{pam}{Adjacent PAM sequence in the genome.}
#'     \item{chr}{Chromosome name.}
#'     \item{pam_site}{Genomic coordinate of the first nucleotide of the PAM.}
#'     \item{strand}{"+" or "-".}
#'     \item{n_mismatches}{Number of mismatches between spacer and protospacer.}
#'     \item{canonical}{Logical. Whether the PAM is canonical for this nuclease.}
#'   }
#'   Returns a zero-row data.frame if no alignments are found.
#'
#' @seealso \code{\link{ensure_bowtie_index}} for index building,
#'          \code{\link{normalise_bowtie_hits}} for schema conversion.
#'
#' @keywords internal
run_bowtie_search <- function(unique_spacers,
                              genome,
                              nuclease = c("Cas9", "Cas12a", "enCas12a"),
                              bowtie_index,
                              max_mismatches = 3L,
                              all_alignments = TRUE,
                              n_max_alignments = 1000L,
                              canonical_only = TRUE,
                              quiet = FALSE) {

  nuclease <- match.arg(nuclease)

  # ---- 1. Validate dependencies ----
  if (!requireNamespace("crisprBowtie", quietly = TRUE)) {
    stop("Package 'crisprBowtie' is required for the Bowtie off-target backend.\n",
         "Install with: BiocManager::install('crisprBowtie')")
  }

  if (!requireNamespace("Rbowtie", quietly = TRUE)) {
    stop("Package 'Rbowtie' is required for the Bowtie off-target backend.\n",
         "Install with: BiocManager::install('Rbowtie')")
  }

  # ---- 2. Validate inputs ----
  if (!is.character(unique_spacers) || length(unique_spacers) == 0) {
    stop("'unique_spacers' must be a non-empty character vector.")
  }

  spacer_lengths <- unique(nchar(unique_spacers))
  if (length(spacer_lengths) != 1) {
    stop("All spacer sequences must be the same length. ",
         "Found lengths: ", paste(sort(spacer_lengths), collapse = ", "))
  }

  if (!inherits(genome, "BSgenome")) {
    stop("'genome' must be a BSgenome object.")
  }

  if (!is.character(bowtie_index) || !nzchar(bowtie_index)) {
    stop("'bowtie_index' must be a non-empty path to a Bowtie index prefix.")
  }

  # Verify index actually exists
  if (!bowtie_index_exists(bowtie_index)) {
    stop("Bowtie index not found at prefix: ", bowtie_index, "\n",
         "Run ensure_bowtie_index() first, or set ot_backend = 'crispritz'.")
  }

  # ---- 3. Clamp mismatches to Bowtie limit ----
  max_mismatches <- as.integer(max_mismatches)
  if (max_mismatches > 3L) {
    warning("Bowtie supports a maximum of 3 mismatches. ",
            "Clamping max_mismatches from ", max_mismatches, " to 3.")
    max_mismatches <- 3L
  }
  if (max_mismatches < 0L) {
    stop("'max_mismatches' must be >= 0.")
  }

  # ---- 4. Get CrisprNuclease object ----
  nuc_obj <- get_crisprbase_nuclease(nuclease)

  # Check spacer length against nuclease expectation
  expected_len <- crisprBase::spacerLength(nuc_obj)
  force_length <- (spacer_lengths != expected_len)

  if (force_length && !quiet) {
    message("Note: spacer length (", spacer_lengths,
            "nt) differs from ", nuclease, " default (",
            expected_len, "nt). Overriding nuclease spacer length.")
  }

  # ---- 5. Run crisprBowtie ----
  if (!quiet) {
    message("Bowtie search: ", length(unique_spacers), " spacers, ",
            max_mismatches, " mismatches, ",
            if (canonical_only) "canonical PAMs only" else "all PAMs")
  }

  t_start <- Sys.time()

  hits_df <- tryCatch({
    crisprBowtie::runCrisprBowtie(
      spacers              = unique_spacers,
      mode                 = "spacer",
      bowtie_index         = bowtie_index,
      bsgenome             = genome,
      crisprNuclease       = nuc_obj,
      canonical            = canonical_only,
      ignore_pam           = FALSE,
      n_mismatches         = max_mismatches,
      all_alignments       = all_alignments,
      n_max_alignments     = as.integer(n_max_alignments),
      force_spacer_length  = force_length,
      verbose              = FALSE # Prevent runCRISPRBowtie verbosity from contaminating pipeline output
    )
  }, error = function(e) {
    stop("crisprBowtie::runCrisprBowtie() failed: ", e$message, "\n",
         "Bowtie index prefix: ", bowtie_index, "\n",
         "Nuclease: ", nuclease, " (", class(nuc_obj)[1], ")\n",
         "Spacer length: ", spacer_lengths, "nt, Mismatches: ", max_mismatches)
  })

  t_elapsed <- difftime(Sys.time(), t_start, units = "secs")

  # ---- 6. Handle NULL/empty return ----
  if (is.null(hits_df) || nrow(hits_df) == 0) {
    if (!quiet) message("Bowtie search completed in ",
                        round(as.numeric(t_elapsed), 1),
                        "s — no alignments found.")
    # Return empty data.frame with expected columns
    return(data.frame(
      spacer       = character(0),
      protospacer  = character(0),
      pam          = character(0),
      chr          = character(0),
      pam_site     = integer(0),
      strand       = character(0),
      n_mismatches = integer(0),
      canonical    = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  if (!quiet) {
    message("Bowtie search completed in ",
            round(as.numeric(t_elapsed), 1), "s — ",
            nrow(hits_df), " alignments found",
            " (", length(unique(hits_df$spacer)), " spacers with hits)")
  }

  return(hits_df)
}

#' Normalise crisprBowtie output to the standardised hit_df schema
#'
#' Converts the data.frame returned by \code{\link{run_bowtie_search}}
#' (i.e. \code{crisprBowtie::runCrisprBowtie()}) into the same column
#' schema produced by \code{\link{parse_crispritz_output}}, so that
#' \code{\link{score_and_aggregate}} can process hits from either backend
#' identically.
#'
#' @param bowtie_df data.frame. Raw output from \code{run_bowtie_search()}.
#'        Expected columns: spacer, protospacer, pam, chr, pam_site, strand,
#'        n_mismatches, canonical.
#' @param guide_map Named list. Maps unique protospacer sequences (names) to
#'        integer vectors of GRanges indices (values), as returned by
#'        \code{\link{deduplicate_protospacers}}.
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a". Used to
#'        determine PAM geometry for protospacer position computation.
#'
#' @return A data.frame with the standardised hit schema:
#'   \describe{
#'     \item{guide_seq}{Character. Bare protospacer sequence (uppercase, no PAM).}
#'     \item{grna_indices}{List-column. Integer vectors of GRanges indices.}
#'     \item{chr}{Character. Chromosome name.}
#'     \item{pos}{Integer. 1-based genomic start of the protospacer.}
#'     \item{strand}{Character. "+" or "-".}
#'     \item{offtarget_seq}{Character. Genomic off-target protospacer (no PAM).}
#'     \item{pam_gen}{Character. PAM sequence found in genome.}
#'     \item{n_mismatches}{Integer. Number of mismatches.}
#'     \item{bulge_type}{Character. Always "X" (Bowtie = mismatch-only).}
#'     \item{bulge_size}{Integer. Always 0.}
#'   }
#'   Returns a zero-row data.frame with correct column types if input is empty.
#'
#' @details
#' \strong{Position derivation.}
#' crisprBowtie reports \code{pam_site} as the leftmost genomic coordinate
#' of the PAM sequence (regardless of strand). The protospacer start is
#' computed from PAM geometry:
#'
#' \itemize{
#'   \item Cas9 (3' PAM), + strand: \code{pam_site - spacer_length}
#'   \item Cas9 (3' PAM), - strand: \code{pam_site + pam_length}
#'   \item Cas12a/enCas12a (5' PAM), + strand: \code{pam_site + pam_length}
#'   \item Cas12a/enCas12a (5' PAM), - strand: \code{pam_site - spacer_length}
#' }
#'
#' This convention should be validated against crisprBowtie test output.
#' The ±10 bp on-target exclusion window in \code{score_and_aggregate}
#' provides tolerance for minor convention differences.
#'
#' @seealso \code{\link{run_bowtie_search}}, \code{\link{parse_crispritz_output}},
#'          \code{\link{score_and_aggregate}}
#'
#' @keywords internal
normalise_bowtie_hits <- function(bowtie_df, guide_map, nuclease) {

  nuclease <- match.arg(nuclease, c("Cas9", "Cas12a", "enCas12a"))

  # ---- 1. Canonical empty schema (matches parse_crispritz_output) ----
  empty_result <- data.frame(
    guide_seq      = character(0),
    grna_indices   = I(list()),
    chr            = character(0),
    pos            = integer(0),
    strand         = character(0),
    offtarget_seq  = character(0),
    pam_gen        = character(0),
    n_mismatches   = integer(0),
    bulge_type     = character(0),
    bulge_size     = integer(0),
    stringsAsFactors = FALSE
  )

  # ---- 2. Handle empty input ----
  if (is.null(bowtie_df) || nrow(bowtie_df) == 0) {
    return(empty_result)
  }

  # ---- 3. Validate expected columns ----
  required_cols <- c("spacer", "protospacer", "pam", "chr",
                     "pam_site", "strand", "n_mismatches")
  missing_cols <- setdiff(required_cols, names(bowtie_df))
  if (length(missing_cols) > 0) {
    stop("crisprBowtie output is missing expected columns: ",
         paste(missing_cols, collapse = ", "),
         "\nAvailable columns: ", paste(names(bowtie_df), collapse = ", "))
  }

  # ---- 4. Extract and clean core columns ----
  guide_seq     <- toupper(as.character(bowtie_df$spacer))
  offtarget_seq <- toupper(as.character(bowtie_df$protospacer))
  pam_gen       <- toupper(as.character(bowtie_df$pam))
  chr           <- as.character(bowtie_df$chr)
  strand        <- as.character(bowtie_df$strand)
  pam_site      <- as.integer(bowtie_df$pam_site)
  n_mismatches  <- as.integer(bowtie_df$n_mismatches)

  # ---- 5. Compute protospacer start from pam_site ----
  # crisprBowtie pam_site = leftmost genomic coordinate of the PAM.
  # Protospacer position depends on PAM side (3' vs 5') and strand.

  spacer_len <- nchar(offtarget_seq)
  pam_len    <- nchar(pam_gen)

  pam_side <- if (nuclease == "Cas9") "3prime" else "5prime"

  if (pam_side == "3prime") {
    # Cas9: PAM is 3' of protospacer on target strand
    #   + strand: [proto] [PAM] → proto_start = pam_site - spacer_len
    #   - strand: [PAM_rc] [proto_rc] → proto_start = pam_site + pam_len
    pos <- ifelse(strand == "+",
                  pam_site - spacer_len,
                  pam_site + pam_len)
  } else {
    # Cas12a/enCas12a: PAM is 5' of protospacer on target strand
    #   + strand: [PAM] [proto] → proto_start = pam_site + pam_len
    #   - strand: [proto_rc] [PAM_rc] → proto_start = pam_site - spacer_len
    pos <- ifelse(strand == "+",
                  pam_site + pam_len,
                  pam_site - spacer_len)
  }

  pos <- as.integer(pos)

  # ---- 6. Map guide sequences to GRanges indices via guide_map ----
  grna_indices <- lapply(guide_seq, function(seq) {
    idx <- guide_map[[seq]]
    if (is.null(idx)) integer(0) else idx
  })

  # Flag unmapped guides (defensive — shouldn't happen with consistent dedup)
  unmapped <- vapply(grna_indices, function(x) length(x) == 0, logical(1))
  if (any(unmapped)) {
    n_unmapped <- sum(unmapped)
    warning(n_unmapped, " Bowtie hit(s) could not be mapped back to input gRNAs. ",
            "These will be excluded from scoring.\n",
            "First unmapped spacer: ", guide_seq[which(unmapped)[1]])

    keep <- !unmapped
    guide_seq     <- guide_seq[keep]
    grna_indices  <- grna_indices[keep]
    chr           <- chr[keep]
    pos           <- pos[keep]
    strand        <- strand[keep]
    offtarget_seq <- offtarget_seq[keep]
    pam_gen       <- pam_gen[keep]
    n_mismatches  <- n_mismatches[keep]
  }

  if (length(guide_seq) == 0) {
    return(empty_result)
  }

  # ---- 7. Assemble standardised hit_df ----
  result <- data.frame(
    guide_seq     = guide_seq,
    grna_indices  = I(grna_indices),
    chr           = chr,
    pos           = pos,
    strand        = strand,
    offtarget_seq = offtarget_seq,
    pam_gen       = pam_gen,
    n_mismatches  = n_mismatches,
    bulge_type    = rep("X", length(guide_seq)),
    bulge_size    = rep(0L, length(guide_seq)),
    stringsAsFactors = FALSE
  )

  return(result)
}
