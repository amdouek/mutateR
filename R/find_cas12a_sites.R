#' Find Cas12a (Cpf1) target sites in exons
#'
#' Scans exon sequences (both strands) for Cas12a gRNA sites (default PAM "TTTV")
#'
#'
#' @param exon_gr GRanges from get_exon_structures(..., output = "GRanges").
#' @param genome  BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg38).
#' @param pam PAM (default "TTTV"). IUPAC ambiguity codes are supported.
#' @param protospacer_length Protospacer length downstream of PAM (default 23 bp).
#'
#' @return GRanges of Cas12a target sites with metadata columns:
#'         exon_rank, protospacer_sequence, pam_sequence, target_sequence,
#'         sequence_context (34 bp), cut_site.
#'
#' @export
find_cas12a_sites <- function(exon_gr,
                              genome,
                              pam = "TTTV",
                              protospacer_length = 23) {

  stopifnot(inherits(exon_gr, "GRanges"),
            "rank" %in% names(mcols(exon_gr)),
            inherits(genome, "BSgenome"))

  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Biostrings package required.")
  GenomeInfoDb::seqlevelsStyle(exon_gr) <- GenomeInfoDb::seqlevelsStyle(genome)

  results <- list()
  pam_len <- nchar(pam)
  pam_fwd <- Biostrings::DNAString(pam)
  pam_rev <- Biostrings::reverseComplement(pam_fwd)

  for (i in seq_along(exon_gr)) {
    gr <- exon_gr[i]
    exon_seq <- Biostrings::getSeq(genome, gr)[[1]]
    exon_len <- nchar(exon_seq)

    ## ----- Forward strand (PAM 5'-TTTV-3') -----
    hits_fwd <- Biostrings::matchPattern(pam_fwd, exon_seq, fixed = FALSE)
    if (length(hits_fwd)) {
      for (h in seq_along(hits_fwd)) {
        pam_start <- start(hits_fwd)[h]
        pam_end   <- end(hits_fwd)[h]

        # ------ Protospacer downstream of PAM ------
        protospacer_start <- pam_end + 1
        protospacer_end   <- pam_end + protospacer_length
        if (protospacer_end > exon_len) next

        pam_seq      <- as.character(Biostrings::subseq(exon_seq, pam_start, pam_end))
        protospacer  <- as.character(Biostrings::subseq(exon_seq, protospacer_start, protospacer_end))
        target       <- paste0(pam_seq, protospacer)

        ## ------- DeepCpf1 34-bp sequence context = 4 nt upstream + 29 downstream from PAM -------
        ctx_start <- max(1, pam_start - 4)
        ctx_end   <- min(exon_len, pam_start + 29)
        seq_context <- as.character(Biostrings::subseq(exon_seq, ctx_start, ctx_end))

        ## ------- Genomic coordinates (relative to exon) -------
        start_genomic <- start(gr) + pam_start - 1
        end_genomic   <- start(gr) + protospacer_end - 1
        cut_site      <- start_genomic + pam_len + 18  # Cas12a cut ~18 bp downstream of PAM

        results[[length(results) + 1]] <- GenomicRanges::GRanges(
          seqnames = seqnames(gr),
          ranges = IRanges::IRanges(start = start_genomic, end = end_genomic),
          strand = "+",
          exon_rank = gr$rank,
          protospacer_sequence = protospacer,
          pam_sequence = pam_seq,
          target_sequence = target,
          cut_site = cut_site,
          sequence_context = seq_context
        )
      }
    }

    ## ----- Reverse strand (PAM revcomp 3'-AAAV-5') -----
    exon_seq_rc <- Biostrings::reverseComplement(exon_seq)
    hits_rev <- Biostrings::matchPattern(pam_rev, exon_seq_rc, fixed = FALSE)
    if (length(hits_rev)) {
      for (h in seq_along(hits_rev)) {
        pam_start <- start(hits_rev)[h]
        pam_end   <- end(hits_rev)[h]

        protospacer_start <- pam_end + 1
        protospacer_end   <- pam_end + protospacer_length
        if (protospacer_end > nchar(exon_seq_rc)) next

        pam_seq      <- as.character(Biostrings::subseq(exon_seq_rc, pam_start, pam_end))
        protospacer  <- as.character(Biostrings::subseq(exon_seq_rc, protospacer_start, protospacer_end))
        target       <- paste0(pam_seq, protospacer)

        ## ------- DeepCpf1 context (on reverse strand): 4 bp upstream and 29 downstream of PAM -------
        ctx_start <- max(1, pam_start - 29)
        ctx_end   <- min(nchar(exon_seq_rc), pam_end + 4)
        raw_ctx   <- Biostrings::subseq(exon_seq_rc, ctx_start, ctx_end)
        seq_context <- as.character(Biostrings::reverseComplement(raw_ctx))

        ## ------- Map to genomic coordinates -------
        pam_genomic_end   <- end(gr) - (pam_start - 1)
        pam_genomic_start <- pam_genomic_end - (pam_len + protospacer_length - 1)
        start_genomic <- pam_genomic_start
        end_genomic   <- pam_genomic_end
        cut_site      <- end_genomic - pam_len - 18

        results[[length(results) + 1]] <- GenomicRanges::GRanges(
          seqnames = seqnames(gr),
          ranges = IRanges::IRanges(start = start_genomic, end = end_genomic),
          strand = "-",
          exon_rank = gr$rank,
          protospacer_sequence = protospacer,
          pam_sequence = pam_seq,
          target_sequence = target,
          cut_site = cut_site,
          sequence_context = seq_context
        )
      }
    }
  }

  if (!length(results)) {
    warning("No Cas12a sites found for PAM ", pam)
    return(NULL)
  }
  do.call(c, results)
}
