#' Find Cas12a (Cpf1) target sites in exons
#'
#' Scans exon sequences (both strands) for Cas12a gRNA sites (default PAM "TTTV")
#' Retrieves full genomic context needed for DeepCpf1 scoring (34 bp) even if sites are near exon boundaries.
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
            inherits(genome, "BSgenome"))

  if (!requireNamespace("Biostrings", quietly = TRUE)) stop("Biostrings required.")
  GenomeInfoDb::seqlevelsStyle(exon_gr) <- GenomeInfoDb::seqlevelsStyle(genome)

  results <- list()
  pam_len <- nchar(pam)
  pam_fwd <- Biostrings::DNAString(pam)
  pam_rev <- Biostrings::reverseComplement(pam_fwd)

  for (i in seq_along(exon_gr)) {
    gr <- exon_gr[i]
    exon_seq <- Biostrings::getSeq(genome, gr)[[1]]
    exon_len <- nchar(exon_seq)

    ## ----- Forward strand -----
    hits_fwd <- Biostrings::matchPattern(pam_fwd, exon_seq, fixed = FALSE)
    if (length(hits_fwd)) {
      for (h in seq_along(hits_fwd)) {
        pam_start <- start(hits_fwd)[h]
        pam_end   <- end(hits_fwd)[h]

        proto_start <- pam_end + 1
        proto_end   <- pam_end + protospacer_length
        if (proto_end > exon_len) next

        # Calculate Genomic Coordinates
        gen_pam_start <- start(gr) + pam_start - 1
        gen_proto_end <- start(gr) + proto_end - 1

        # Calculate Context Coordinates (34 bp window)
        # 4 bp upstream of PAM start
        ctx_gen_start <- gen_pam_start - 4
        # End = Start + 33 (Total 34 bp)
        ctx_gen_end   <- ctx_gen_start + 33

        # Safe fetch from genome (handles intron/flank)
        seq_context <- tryCatch({
          as.character(Biostrings::getSeq(genome,
                                          GenomicRanges::GRanges(seqnames(gr),
                                                                 IRanges(ctx_gen_start, ctx_gen_end),
                                                                 strand="+")))
        }, error = function(e) NA_character_)

        # Skip if fetch failed (e.g. chromosome edge)
        if (is.na(seq_context) || nchar(seq_context) != 34) next

        pam_seq      <- as.character(Biostrings::subseq(exon_seq, pam_start, pam_end))
        protospacer  <- as.character(Biostrings::subseq(exon_seq, proto_start, proto_end))

        results[[length(results) + 1]] <- GenomicRanges::GRanges(
          seqnames = seqnames(gr),
          ranges = IRanges::IRanges(start = gen_pam_start, end = gen_proto_end),
          strand = "+",
          exon_rank = gr$rank,
          protospacer_sequence = protospacer,
          pam_sequence = pam_seq,
          target_sequence = paste0(pam_seq, protospacer),
          cut_site = gen_pam_start + pam_len + 18,
          sequence_context = seq_context
        )
      }
    }

    ## ----- Reverse strand -----
    exon_seq_rc <- Biostrings::reverseComplement(exon_seq)
    hits_rev <- Biostrings::matchPattern(pam_rev, exon_seq_rc, fixed = FALSE)
    if (length(hits_rev)) {
      for (h in seq_along(hits_rev)) {
        pam_start_rc <- start(hits_rev)[h]
        pam_end_rc   <- end(hits_rev)[h]

        proto_start_rc <- pam_end_rc + 1
        proto_end_rc   <- pam_end_rc + protospacer_length
        if (proto_end_rc > nchar(exon_seq_rc)) next

        # Genomic Mapping (Reverse Strand)
        # RC index 1 = Genomic End. RC Index k = Genomic End - (k-1)
        gen_pam_start_5p <- end(gr) - (pam_start_rc - 1) # 5' of PAM on minus strand (high coord)
        gen_proto_end_3p <- end(gr) - (proto_end_rc - 1) # 3' of Proto on minus strand (low coord)

        # Context on Minus Strand
        # We need 4bp upstream (higher coord), 29bp downstream (lower coord) relative to PAM 5'
        # GRanges handles the strand inversion. We just need to define the span.
        # High Coord (Genomic Start of Context on Minus) = PAM_5p + 4
        # Low Coord (Genomic End of Context on Minus) = High - 33

        ctx_gen_high <- gen_pam_start_5p + 4
        ctx_gen_low  <- ctx_gen_high - 33

        seq_context <- tryCatch({
          as.character(Biostrings::getSeq(genome,
                                          GenomicRanges::GRanges(seqnames(gr),
                                                                 IRanges(ctx_gen_low, ctx_gen_high),
                                                                 strand="-")))
        }, error = function(e) NA_character_)

        if (is.na(seq_context) || nchar(seq_context) != 34) next

        pam_seq      <- as.character(Biostrings::subseq(exon_seq_rc, pam_start_rc, pam_end_rc))
        protospacer  <- as.character(Biostrings::subseq(exon_seq_rc, proto_start_rc, proto_end_rc))

        results[[length(results) + 1]] <- GenomicRanges::GRanges(
          seqnames = seqnames(gr),
          ranges = IRanges::IRanges(start = gen_proto_end_3p, end = gen_pam_start_5p),
          strand = "-",
          exon_rank = gr$rank,
          protospacer_sequence = protospacer,
          pam_sequence = pam_seq,
          target_sequence = paste0(pam_seq, protospacer),
          cut_site = gen_pam_start_5p - pam_len - 18,
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
