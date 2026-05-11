#' @title Find Cas9 target sites in exons
#'
#' @description Scans exon sequences on both strands for NGG PAMs.
#' Excludes gRNA sequences where the flanking -4/+3 nucleotides would run out of an exon boundary (for on-target scoring purposes).
#' Outputs both + and - strand hits, each with a correctly oriented 30‑nt
#' sequence_context = [4-nt flank][20-nt protospacer][3-nt PAM][3-nt flank].
#'
#' @import Biostrings
#'
#'
#' @param exon_gr GRanges from get_exon_structures(..., output="GRanges").
#'        Must contain metadata column `rank`.
#' @param genome BSgenome object (e.g. BSgenome.Hsapiens.UCSC.hg38).
#' @param pam PAM (default "NGG"). IUPAC ambiguity codes are supported.
#' @param protospacer_length Protospacer length upstream of PAM (default 20 bp for Cas9).
#'
#' @return GRanges object of candidate scorable Cas9 binding sites with relevant metadata.
#'
#' @export

find_cas9_sites <- function(exon_gr,
                            genome,
                            pam = "NGG",
                            protospacer_length = 20) {

  stopifnot(inherits(exon_gr, "GRanges"),
            "rank" %in% names(mcols(exon_gr)),
            inherits(genome, "BSgenome"))

  GenomeInfoDb::seqlevelsStyle(exon_gr) <- GenomeInfoDb::seqlevelsStyle(genome)
  pam_seq <- Biostrings::DNAString(pam)
  hits_all <- list()

  # ----- Helper function: Safe substring extraction -----
  safe_subseq <- function(x, start, end)
    Biostrings::subseq(x,
                       start = max(start, 1),
                       end   = min(end, nchar(x)))

  for (i in seq_along(exon_gr)) {
    gr <- exon_gr[i]
    exon_seq <- Biostrings::getSeq(genome, gr)[[1]]
    exon_len <- nchar(exon_seq)

    for (strand_sign in c("+", "-")) {
      if (strand_sign == "+") {
        seq_to_scan <- exon_seq
        offset_func <- function(hit_start, hit_end)
          list(start    = start(gr) + hit_start - protospacer_length - 1,
               end      = start(gr) + hit_end - 1,
               cut_site = start(gr) + hit_start - protospacer_length - 1 +
                            protospacer_length - 3)
        #  Scan the forward strand (of the transcript, not genomic).
        #  On the + strand the protospacer occupies positions
        #  [hit_start - proto_len, hit_start - 1] of exon_seq; the genomic
        #  range of the (proto + PAM) hit therefore runs from
        #  start(gr) + hit_start - proto_len - 1 (proto 5') to
        #  start(gr) + hit_end - 1 (PAM 3'). The blunt SpCas9 cut sits 3 nt
        #  upstream of the PAM, i.e. at proto position 18 from the 5' end.
      } else {
        seq_to_scan <- Biostrings::reverseComplement(exon_seq)
        offset_func <- function(hit_start, hit_end) {
          # `pam_end_genomic` is named for symmetry with the + strand block
          # but in fact represents the *high* genomic coord of the PAM on the
          # minus strand — i.e. the PAM 5' nucleotide when reading the target
          # strand 5'→3'. On the minus strand target reads high→low in
          # genomic coords, so the (proto + PAM) range lays out as:
          #   genomic low  = PAM 3' on target = pam_end_genomic - pam_len + 1
          #   genomic high = proto 5' on target = pam_end_genomic + proto_len
          # (Previous versions placed the range pam_len + proto_len - 1 bp
          # *below* the PAM, which put it 22 bp downstream of the actual
          # protospacer; cut_site was also off by 8 bp on the minus strand.)
          pam_end_genomic   <- end(gr) - (hit_start - 1)
          pam_len_n         <- nchar(pam)
          range_lo_genomic  <- pam_end_genomic - pam_len_n + 1L
          range_hi_genomic  <- pam_end_genomic + protospacer_length
          # Blunt cut 3 nt 5' of the PAM on target = 3 nt above the PAM in
          # genomic coords on the minus strand.
          cut_genomic       <- pam_end_genomic + 3L
          list(start    = range_lo_genomic,
               end      = range_hi_genomic,
               cut_site = cut_genomic)
        }
      }

      hits <- Biostrings::matchPattern(pam_seq, seq_to_scan, fixed = FALSE)
      if (!length(hits)) next

      for (h in seq_along(hits)) {
        pam_start <- Biostrings::start(hits)[h]
        pam_end   <- Biostrings::end(hits)[h]

        # Ensure complete protospacer and flanks fit within exon sequence (for on-target scoring purposes) <- To consider: Make this a user-definable behaviour in case the user is happy with protospacers extending into introns. However, would require separate retrieval of intronic sequences which is currently not supported.
        left_bound  <- pam_start - protospacer_length - 4     # inclusive
        right_bound <- pam_end + 3                            # inclusive
        if (left_bound < 1 || right_bound > exon_len) next     # skip out-of-exon-bound cases

        # ----- Extract protospacer and PAM in this scanning orientation -----
        protospacer <- Biostrings::subseq(seq_to_scan,
                                          pam_start - protospacer_length,
                                          pam_start - 1)
        pam_seq_str <- Biostrings::subseq(seq_to_scan, pam_start, pam_end)
        # flanks
        left_flank  <- Biostrings::subseq(seq_to_scan,
                                          pam_start - protospacer_length - 4,
                                          pam_start - protospacer_length - 1)
        right_flank <- Biostrings::subseq(seq_to_scan, pam_end + 1, pam_end + 3)

        ctx <- paste0(left_flank,
                      protospacer,
                      pam_seq_str,
                      right_flank) #  Assembles full gRNA sequence context

        # ----- Genomic coordinate mapping -----
        coords <- offset_func(pam_start, pam_end)

        hit <- GRanges(
          seqnames = seqnames(gr),
          IRanges(start = coords$start, end = coords$end),
          strand = strand_sign,
          exon_rank = gr$rank,
          protospacer_sequence = as.character(protospacer),
          pam_sequence          = as.character(pam_seq_str),
          target_sequence       = as.character(paste0(protospacer, pam_seq_str)),
          cut_site              = coords$cut_site,
          sequence_context      = as.character(ctx)
        )

        hits_all[[length(hits_all) + 1]] <- hit
      }
    }
  }

  if (!length(hits_all)) {
    warning("No Cas9 sites found for PAM ", pam)
    return(NULL)
  }
  do.call(c, hits_all)
}
