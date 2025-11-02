#' Assemble and annotate gRNA pairs for exon‑flanking deletions
#'
#' For each *phase‑compatible* exon pair (E5′, E3′),
#' gRNAs are selected from the immediately flanking exons
#' (E5′ + 1 and E3′ – 1).
#'
#' The output merges each gRNA pair’s protospacer/PAM sequences with the
#' corresponding on‑target scores retrieved directly from the nested
#' \code{ontarget_score} data.frames in the GRanges metadata.
#'
#' Each 30‑nt scoring sequence is internally trimmed by removing the
#' first 4 and last 6 bases (3‑nt PAM + 3‑nt flank) to recover the
#' 20‑nt protospacer used for matching.
#'
#' @param grna_gr        GRanges returned by \link{filter_valid_grnas},
#'                       containing metadata column \code{ontarget_score}
#'                       (a data.frame with columns \code{sequence} and \code{score}).
#' @param exon_gr        GRanges from \link{get_exon_structures}(output="GRanges").
#' @param transcript_id  Character. Ensembl transcript ID.
#' @param species        Character. Ensembl species short code (e.g. "hsapiens").
#' @param score_cutoff   Numeric. Minimum acceptable on‑target score for each gRNA (default 0.6).
#'
#' @return A \code{data.frame} where each row represents a candidate
#'         exon‑flanking gRNA pair with protospacer/PAM sequences,
#'         matched numeric on‑target scores, frameshift/PTC metadata,
#'         domain annotation, and a logical \code{recommended} flag.
#'
#' @seealso \link{filter_valid_grnas}, \link{check_exon_phase},
#'          \link{check_frameshift_ptc}, \link{map_protein_domains}
#'
#' @examples
#' \dontrun{
#' tp53_pairs_df <- assemble_grna_pairs(
#'   grna_gr       = tp53_valid_grnas,
#'   exon_gr       = tp53_exons_gr,
#'   transcript_id = canonical_tx,
#'   species       = "hsapiens",
#'   score_cutoff  = 0.6
#' )
#' head(tp53_pairs_df)
#' }
#' @export
assemble_grna_pairs <- function(grna_gr,
                                exon_gr,
                                transcript_id,
                                species,
                                score_cutoff = 0.5) {

  stopifnot(inherits(grna_gr, "GRanges"))
  message("Assembling gRNA pairs for exon‑flanking deletions...")

  ## ---- Handling for single‑ or two‑exon transcripts --------------------
  n_exons <- length(exon_gr)
  if (n_exons <= 2) {
    message("Single-exon/two-exon gene detected: constructing intragenic deletion pairs.")

    # ---- Flatten ontarget_score safely to numeric -------------------
    raw_scores <- mcols(grna_gr)$ontarget_score

    if (is.numeric(raw_scores)) {
      numeric_scores <- raw_scores

    } else if (is.data.frame(raw_scores) && "score" %in% names(raw_scores)) {
      numeric_scores <- suppressWarnings(as.numeric(raw_scores$score))

    } else if (inherits(raw_scores, "DataFrame") && "score" %in% names(raw_scores)) {
      numeric_scores <- as.numeric(raw_scores$score)

    } else if (inherits(raw_scores, "DataFrame")) {
      # If it's a one‑column DataFrame named something else
      numeric_scores <- suppressWarnings(as.numeric(raw_scores[[1]]))

    } else if (is.list(raw_scores)) {
      extract_num <- function(x) {
        if (is.numeric(x)) return(x[1])
        if (is.data.frame(x) && "score" %in% names(x)) return(as.numeric(x$score[1]))
        if (is.data.frame(x)) return(suppressWarnings(as.numeric(x[[1]][1])))
        if (inherits(x, "DataFrame")) return(suppressWarnings(as.numeric(x[[1]][1])))
        if (is.list(x) && length(x) == 1) return(Recall(x[[1]]))
        suppressWarnings(as.numeric(x[1]))
      }
      numeric_scores <- vapply(raw_scores, extract_num, numeric(1))

    } else {
      numeric_scores <- suppressWarnings(as.numeric(raw_scores))
    }

    numeric_scores[is.nan(numeric_scores)] <- NA_real_

    # Build tabular metadata
    site_df <- as.data.frame(grna_gr)
    # derive protospacer if absent
    if (!"protospacer_sequence" %in% names(mcols(grna_gr))) {
      if ("target_sequence" %in% names(mcols(grna_gr))) {
        site_df$protospacer_sequence <- mcols(grna_gr)$target_sequence
      } else if ("sequence_context" %in% names(mcols(grna_gr))) {
        sc <- as.character(mcols(grna_gr)$sequence_context)
        site_df$protospacer_sequence <- substring(sc, 5, nchar(sc) - 6)
      } else {
        stop("Cannot derive protospacer sequences from provided GRanges object.")
      }
    } else {
      site_df$protospacer_sequence <- mcols(grna_gr)$protospacer_sequence
    }

    # exon rank assignment if absent
    if (!"exon_rank" %in% names(mcols(grna_gr))) {
      site_df$exon_rank <- 1L
    } else {
      site_df$exon_rank <- mcols(grna_gr)$exon_rank
    }

    site_df$pos <- GenomicRanges::start(grna_gr)
    site_df$score <- numeric_scores

    if (nrow(site_df) < 2) {
      warning("Fewer than two scored gRNAs available for pairing.")
      return(NULL)
    }

    # build all guide pairs
    comb_idx <- utils::combn(seq_len(nrow(site_df)), 2)
    intragenic <- data.frame(
      protospacer_sequence_5p = site_df$protospacer_sequence[comb_idx[1,]],
      protospacer_sequence_3p = site_df$protospacer_sequence[comb_idx[2,]],
      ontarget_score_5p       = site_df$score[comb_idx[1,]],
      ontarget_score_3p       = site_df$score[comb_idx[2,]],
      exon_5p                 = site_df$exon_rank[comb_idx[1,]],
      exon_3p                 = site_df$exon_rank[comb_idx[2,]],
      del_size                = abs(site_df$pos[comb_idx[2,]] - site_df$pos[comb_idx[1,]]),
      transcript_id           = transcript_id,
      stringsAsFactors = FALSE
    )

    # rank and recommend
    intragenic <- intragenic[order(-intragenic$del_size), ]
    intragenic$recommended <-
      with(intragenic,
           ontarget_score_5p >= score_cutoff &
             ontarget_score_3p >= score_cutoff)

    message("Generated ", nrow(intragenic), " intragenic deletion pairs; ",
            sum(intragenic$recommended), " meet score cutoff.")

    # Keep only recommended pairs for final output
    intragenic <- intragenic[intragenic$recommended == TRUE, ]
    message("Returning ", nrow(intragenic), " recommended intragenic pairs.")

    # Wrap in named list to mimic standard output structure
    return(list(
      pairs = intragenic,
      intragenic_mode = TRUE
    ))
  }

  ## ----- Default multi-exon logic -----
  # if intragenic_mode flag exists, stop further assembly logic
  if (exists("intragenic_mode", inherits = FALSE)) return(pairs)

  exon_meta <- as.data.frame(mcols(exon_gr))
  exon_meta$rank <- seq_len(nrow(exon_meta))
  comp_pairs <- check_exon_phase(exon_meta, include_contiguous = FALSE)
  comp_pairs <- subset(comp_pairs, compatible == TRUE)
  if (nrow(comp_pairs) == 0) {
    warning("No phase‑compatible exon pairs found.")
    return(NULL)
  }

  # Assess frameshift/PTC status
  fs_list <- lapply(seq_len(nrow(comp_pairs)), function(i)
    with(comp_pairs[i, ],
         check_frameshift_ptc(exon_meta,
                              exon_5p = exon_5p,
                              exon_3p = exon_3p)))
  fs_df <- do.call(rbind, lapply(fs_list, as.data.frame))
  pair_info <- cbind(comp_pairs, fs_df)

  # Keep phase-compatible or terminal‑exon cases
  pair_info <- subset(pair_info, compatible | terminal_exon_case)

  # Define flanking exons
  n_exons <- nrow(exon_meta)
  pair_info$target_5p <- pair_info$exon_5p + 1
  pair_info$target_3p <- pair_info$exon_3p - 1
  pair_info <- subset(pair_info,
                      target_5p > 0 & target_3p <= n_exons & target_5p < target_3p)
  if (nrow(pair_info) == 0) {
    warning("No eligible flanking exon pairs after bounds filtering.")
    return(NULL)
  }


  ## ---- Build flat on‑target‑score lookup -----------------------------
  # ---- Build flat on-target-score lookup ----------------------
  message("Flattening on-target scores from GRanges ...")

  # ensure dplyr is loaded internally
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' needed for dataframe manipulation.")

  ontarget_list <- mcols(grna_gr)$ontarget_score

  if (is.list(ontarget_list) && all(vapply(ontarget_list, is.data.frame, logical(1)))) {
    score_lookup <- dplyr::bind_rows(ontarget_list)
  } else if (inherits(ontarget_list, "DataFrame")) {
    score_lookup <- as.data.frame(ontarget_list)
  } else if (is.data.frame(ontarget_list)) {
    score_lookup <- ontarget_list
  } else {
    score_lookup <- data.frame(
      sequence = as.character(mcols(grna_gr)$sequence_context),
      score = suppressWarnings(as.numeric(mcols(grna_gr)$ontarget_score))
    )
  }

  # Clean and standardise
  score_lookup <- dplyr::mutate(score_lookup,
                                protospacer_sequence = substring(sequence, 5, nchar(sequence) - 6))
  score_lookup <- dplyr::select(score_lookup,
                                protospacer_sequence, score)
  score_lookup <- dplyr::distinct(score_lookup)
  score_lookup$score <- suppressWarnings(as.numeric(score_lookup$score))
  score_lookup <- score_lookup[!is.na(score_lookup$score), , drop = FALSE]


  ## ---- Guide metadata + exon rank mapping ----------------------------
  grna_df <- as.data.frame(grna_gr)
  if (!"exon_rank" %in% names(grna_df)) {
    hits <- GenomicRanges::findOverlaps(grna_gr, exon_gr)
    grna_df$exon_rank <- NA_integer_
    grna_df$exon_rank[queryHits(hits)] <-
      exon_meta$rank[subjectHits(hits)]
  }


  ## ---- Domain annotations (optional - improve when final domain annotation source decided) ---------------------------------
  domain_df <- tryCatch(map_protein_domains(transcript_id, species),
                        error = \(e) NULL)
  if (is.null(domain_df)) domain_df <- data.frame()


  ## ---- Build all exon‑flanking gRNA pair combinations ----------------
  results <- list()
  for (i in seq_len(nrow(pair_info))) {
    e5 <- pair_info$target_5p[i]
    e3 <- pair_info$target_3p[i]
    g5 <- subset(grna_df, exon_rank == e5)
    g3 <- subset(grna_df, exon_rank == e3)
    if (nrow(g5) == 0 || nrow(g3) == 0) next

    comb <- merge(g5, g3, by = NULL, suffixes = c("_5p", "_3p"))
    comb$exon_5p <- e5
    comb$exon_3p <- e3
    comb$upstream_pair <- pair_info$exon_5p[i]
    comb$downstream_pair <- pair_info$exon_3p[i]
    comb$compatible <- pair_info$compatible[i]
    comb$frameshift <- pair_info$frameshift[i]
    comb$ptc_flag <- pair_info$ptc_flag[i]
    comb$terminal_exon_case <- pair_info$terminal_exon_case[i]

    dom_overlap <- NULL
    if (nrow(domain_df) > 0 && "exon_rank" %in% names(domain_df))
      dom_overlap <- subset(domain_df, exon_rank %in% seq(e5, e3))
    comb$domains <- if (!is.null(dom_overlap) && nrow(dom_overlap) > 0)
      paste(unique(dom_overlap$domain_desc), collapse = "; ")
    else NA_character_

    results[[length(results) + 1]] <- comb
  }

  if (!length(results)) {
    warning("No valid gRNA pairs identified after flanking‑exon search.")
    return(NULL)
  }
  out <- do.call(rbind, results)

  # Ensure logical consistency between flags
  out$ptc_flag <- ifelse(out$terminal_exon_case, FALSE, out$ptc_flag)
  out$frameshift <- ifelse(out$terminal_exon_case, TRUE, out$frameshift)

  ## ---- Join numeric scores ------------------------------------------
  out <- merge(out, score_lookup,
               by.x = "protospacer_sequence_5p",
               by.y = "protospacer_sequence",
               all.x = TRUE,
               suffixes = c("", "_5p"))
  names(out)[names(out) == "score_5p"] <- "ontarget_score_5p"

  out <- merge(out, score_lookup,
               by.x = "protospacer_sequence_3p",
               by.y = "protospacer_sequence",
               all.x = TRUE,
               suffixes = c("", "_3p"))
  names(out)[names(out) == "score_3p"] <- "ontarget_score_3p"

  ## ---- 'Recommended' flag (score cutoff still required) ---------------
  out$recommended <- with(
    out,
    suppressWarnings(as.numeric(ontarget_score_5p)) >= score_cutoff &
      suppressWarnings(as.numeric(ontarget_score_3p)) >= score_cutoff
  )

  # Optional helper flag to annotate tolerated but low‑score terminal‑exon cases
  out$terminal_tolerated <- out$terminal_exon_case & !out$recommended

  ## ---- Final tidy output --------------------------------------------
  keep_cols <- c("upstream_pair","downstream_pair","exon_5p","exon_3p",
                 "compatible","frameshift","ptc_flag","terminal_exon_case",
                 "protospacer_sequence_5p","pam_sequence_5p",
                 "ontarget_score_5p","protospacer_sequence_3p",
                 "pam_sequence_3p","ontarget_score_3p",
                 "domains","recommended")
  keep_cols <- intersect(keep_cols, names(out))
  out <- out[, keep_cols, drop = FALSE]

  message("Generated ", nrow(out),
          " candidate exon‑flanking gRNA pairs.")

  return(out)
}
