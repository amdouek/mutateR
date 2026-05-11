#' @title Find gRNA target sites for a custom NucleaseSpec
#'
#' @description Generic site finder driven by a user-supplied
#' \code{\link{nuclease_spec}}. Used by the custom-nuclease pipeline path
#' (\code{\link{run_mutateR_custom}}, and \code{\link{run_mutateR}} when called
#' with a custom spec). Does NOT replace the canonical \code{\link{find_cas9_sites}}
#' and \code{\link{find_cas12a_sites}} functions, which continue to drive the
#' canonical Cas9 / Cas12a / enCas12a pipelines.
#'
#' @details
#' \strong{Coordinate convention.} Every site is computed in
#' "target-strand reading" (5'→3' on the strand the gRNA hybridises to). For
#' minus-strand sites this is the reverse complement of the exon sequence.
#' This unifies the two PAM topologies into a single arithmetic and produces
#' correct genomic ranges regardless of strand.
#'
#' \strong{Cut-site geometry.} The spec carries \code{cut_offset_top} and
#' \code{cut_offset_bottom}, measured relative to the protospacer 3' end on
#' the target strand. For each site the finder emits:
#' \itemize{
#'   \item \code{cut_site_top}: genomic coordinate of the nt immediately 3' of
#'         the top-strand cleavage in target-strand reading.
#'   \item \code{cut_site_bottom}: same for the bottom-strand cleavage.
#'   \item \code{cut_site}: the integer midpoint
#'         \code{floor((cut_site_top + cut_site_bottom) / 2)}. Equal to both
#'         endpoints for blunt cutters.
#' }
#' Downstream pair assembly consumes \code{cut_site}; the two endpoints remain
#' available for users who need the full stagger geometry (e.g. for
#' overhang-aware repair design or visualisation).
#'
#' @param exon_gr GRanges from \code{\link{get_exon_structures}}(output="GRanges").
#'        Must contain metadata column \code{rank}.
#' @param genome BSgenome object (e.g. BSgenome.Hsapiens.UCSC.hg38).
#' @param nuclease_spec A \code{NucleaseSpec} object built with
#'        \code{\link{nuclease_spec}}. Canonical strings (\code{"Cas9"} etc.)
#'        are accepted and resolved internally, but for canonical workflows
#'        prefer the dedicated \code{\link{find_cas9_sites}} /
#'        \code{\link{find_cas12a_sites}} entry points.
#' @param require_full_context Logical or NULL. If \code{TRUE}, sites are
#'        skipped when a complete \code{sequence_context} cannot be assembled
#'        within the exon boundaries. If \code{FALSE} (default for custom
#'        specs), sites at exon boundaries are retained and
#'        \code{sequence_context} is fetched from the genome where possible
#'        (intronic flanks allowed); set to \code{NA_character_} when fetch
#'        fails (e.g. chromosome edge). If \code{NULL} (default), behaviour
#'        depends on the spec: \code{TRUE} for canonical specs (preserves
#'        existing scoring-aware behaviour), \code{FALSE} for custom specs.
#'
#' @return A GRanges of candidate sites with metadata columns:
#'   \describe{
#'     \item{\code{exon_rank}}{Integer.}
#'     \item{\code{protospacer_sequence}, \code{pam_sequence}, \code{target_sequence}}{Character.}
#'     \item{\code{sequence_context}}{Character; total width = \code{proto_len + pam_len + 7}.}
#'     \item{\code{cut_site}}{Integer (midpoint of top/bottom cuts in genomic coords).}
#'     \item{\code{cut_site_top}, \code{cut_site_bottom}}{Integer (genomic coords of
#'           the two strand cuts; equal for blunt cutters).}
#'     \item{\code{nuclease}}{Character (the spec's \code{name}).}
#'   }
#'   Returns NULL with a warning when no sites are found.
#'
#' @seealso \code{\link{nuclease_spec}}, \code{\link{run_mutateR_custom}}.
#'
#' @export
find_custom_grna_sites <- function(exon_gr,
                                   genome,
                                   nuclease_spec,
                                   require_full_context = NULL) {

  stopifnot(inherits(exon_gr, "GRanges"),
            "rank" %in% names(GenomicRanges::mcols(exon_gr)),
            inherits(genome, "BSgenome"))

  spec <- resolve_nuclease(nuclease_spec)

  if (is.null(require_full_context)) {
    require_full_context <- isTRUE(spec$is_canonical)
  }
  if (!is.logical(require_full_context) || length(require_full_context) != 1L) {
    stop("`require_full_context` must be TRUE, FALSE, or NULL.")
  }

  GenomeInfoDb::seqlevelsStyle(exon_gr) <- GenomeInfoDb::seqlevelsStyle(genome)

  pam_len    <- spec$pam_length
  proto_len  <- spec$protospacer_length
  pam_side   <- spec$pam_side
  offset_top <- spec$cut_offset_top
  offset_bot <- spec$cut_offset_bottom
  pam_pat    <- Biostrings::DNAString(spec$pam)

  # Context-window dimensions follow the existing convention:
  #   Cas9-style  (3' PAM): [4-nt flank][proto][PAM][3-nt flank]
  #   Cas12a-style (5' PAM): [4-nt flank][PAM][proto][3-nt flank]
  # Same total width either way; only the layout of PAM vs. proto differs.
  ctx_left  <- 4L
  ctx_right <- 3L
  ctx_total <- proto_len + pam_len + ctx_left + ctx_right

  results <- list()

  for (i in seq_along(exon_gr)) {

    gr       <- exon_gr[i]
    exon_seq <- Biostrings::getSeq(genome, gr)[[1]]
    exon_len <- nchar(exon_seq)

    for (strand_sign in c("+", "-")) {

      seq_to_scan <- if (strand_sign == "+") {
        exon_seq
      } else {
        Biostrings::reverseComplement(exon_seq)
      }

      hits <- Biostrings::matchPattern(pam_pat, seq_to_scan, fixed = FALSE)
      if (!length(hits)) next

      for (h in seq_along(hits)) {

        pam_local_start <- Biostrings::start(hits)[h]
        pam_local_end   <- Biostrings::end(hits)[h]

        # ---- Protospacer local coords (target-strand reading) ----
        if (pam_side == "3prime") {
          proto_local_start <- pam_local_start - proto_len
          proto_local_end   <- pam_local_start - 1L
          hit_local_lo      <- proto_local_start
          hit_local_hi      <- pam_local_end
        } else {  # 5prime
          proto_local_start <- pam_local_end + 1L
          proto_local_end   <- pam_local_end + proto_len
          hit_local_lo      <- pam_local_start
          hit_local_hi      <- proto_local_end
        }

        # Reject if protospacer extends past the exon
        if (proto_local_start < 1L || proto_local_end > exon_len) next

        ctx_local_lo <- hit_local_lo - ctx_left
        ctx_local_hi <- hit_local_hi + ctx_right
        full_context_in_exon <- (ctx_local_lo >= 1L && ctx_local_hi <= exon_len)

        if (require_full_context && !full_context_in_exon) next

        # ---- Sequence extraction (target-strand reading) ----
        protospacer <- as.character(Biostrings::subseq(seq_to_scan,
                                                       proto_local_start,
                                                       proto_local_end))
        pam_seq_str <- as.character(Biostrings::subseq(seq_to_scan,
                                                       pam_local_start,
                                                       pam_local_end))
        target_seq  <- if (pam_side == "3prime") {
          paste0(protospacer, pam_seq_str)
        } else {
          paste0(pam_seq_str, protospacer)
        }

        # ---- Genomic coordinate mapping ----
        # Local target-reading position k -> genomic:
        #   + strand: start(gr) + k - 1
        #   - strand: end(gr) - k + 1
        if (strand_sign == "+") {
          hit_gen_lo <- start(gr) + hit_local_lo - 1L
          hit_gen_hi <- start(gr) + hit_local_hi - 1L
        } else {
          hit_gen_lo <- end(gr) - hit_local_hi + 1L
          hit_gen_hi <- end(gr) - hit_local_lo + 1L
        }

        # ---- Sequence context ----
        if (full_context_in_exon) {
          ctx <- as.character(Biostrings::subseq(seq_to_scan,
                                                 ctx_local_lo,
                                                 ctx_local_hi))
        } else {
          if (strand_sign == "+") {
            ctx_gen_lo <- start(gr) + ctx_local_lo - 1L
            ctx_gen_hi <- start(gr) + ctx_local_hi - 1L
          } else {
            ctx_gen_lo <- end(gr) - ctx_local_hi + 1L
            ctx_gen_hi <- end(gr) - ctx_local_lo + 1L
          }
          ctx <- tryCatch({
            ctx_gr <- GenomicRanges::GRanges(
              seqnames = seqnames(gr),
              ranges   = IRanges::IRanges(start = max(ctx_gen_lo, 1L),
                                          end   = ctx_gen_hi),
              strand   = strand_sign
            )
            ctx_chr <- as.character(Biostrings::getSeq(genome, ctx_gr))
            if (nchar(ctx_chr) != ctx_total) NA_character_ else ctx_chr
          }, error = function(e) NA_character_)
        }

        # ---- Cut-site geometry ----
        # Cut nt in target-reading local coords:
        #   proto_local_end + offset + 1
        # (offset is signed, relative to the protospacer 3' end)
        cut_local_top <- proto_local_end + offset_top + 1L
        cut_local_bot <- proto_local_end + offset_bot + 1L

        if (strand_sign == "+") {
          cut_gen_top <- start(gr) + cut_local_top - 1L
          cut_gen_bot <- start(gr) + cut_local_bot - 1L
        } else {
          cut_gen_top <- end(gr) - cut_local_top + 1L
          cut_gen_bot <- end(gr) - cut_local_bot + 1L
        }
        cut_site_mid <- as.integer(floor((cut_gen_top + cut_gen_bot) / 2))

        results[[length(results) + 1]] <- GenomicRanges::GRanges(
          seqnames = seqnames(gr),
          ranges   = IRanges::IRanges(start = hit_gen_lo, end = hit_gen_hi),
          strand   = strand_sign,
          exon_rank            = gr$rank,
          protospacer_sequence = protospacer,
          pam_sequence         = pam_seq_str,
          target_sequence      = target_seq,
          cut_site             = cut_site_mid,
          cut_site_top         = as.integer(cut_gen_top),
          cut_site_bottom      = as.integer(cut_gen_bot),
          sequence_context     = ctx,
          nuclease             = spec$name
        )
      }
    }
  }

  if (!length(results)) {
    warning("No ", spec$name, " sites found for PAM ", spec$pam)
    return(NULL)
  }

  do.call(c, results)
}
