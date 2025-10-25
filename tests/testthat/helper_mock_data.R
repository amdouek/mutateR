# mock one-exon gene for local testing, uses DNAString instead of BSgenome
mock_seq <- Biostrings::DNAString(paste(rep("ATGCGTACGTAGCTAGCTAGCGTAGCTAGCGATCG", 3), collapse = ""))
mock_gr <- GenomicRanges::GRanges("chr1:100-250:+")
mcols(mock_gr)$rank <- 1L
mcols(mock_gr)$start_phase <- 0L
mcols(mock_gr)$end_phase <- 0L
mcols(mock_gr)$cds_start <- 100L
mcols(mock_gr)$cds_end <- 250L
mcols(mock_gr)$exon_cds_length <- 151L