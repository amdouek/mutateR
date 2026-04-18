#' @title Off-target scoring utility functions
#'
#' @description Shared utility functions used across off-target scoring backends
#' (CRISPRitz, Bowtie, hybrid). Includes platform-safe file I/O and
#' protospacer deduplication logic.
#'
#' @keywords internal
#' @name ot_utils
NULL

#' Write text lines to a file with guaranteed LF (Unix) line endings
#'
#' Drop-in replacement for \code{writeLines()} that avoids CRLF on Windows.
#' All content is written in binary mode with explicit LF separators,
#' ensuring compatibility with Unix-native tools (CRISPRitz, Bowtie) when
#' running under WSL or Linux.
#'
#' NOTE: This function escapes $ characters for shell safety.
#' Do NOT use for files containing literal $ references (e.g. awk scripts).
#' For awk scripts, use \code{write_awk_script()} which writes via a binary
#' connection with no character transformation.
#'
#' @param text Character vector. Lines to write (one element per line).
#' @param path Character. Output file path.
#'
#' @return Invisible \code{path}.
#' @keywords internal
write_unix_lines <- function(text, path) {
  text <- gsub("\r", "", text, fixed = TRUE)
  out <- paste0(paste(text, collapse = "\n"), "\n")
  writeBin(charToRaw(out), path)
  invisible(path)
}


#' Strip carriage-return bytes from an existing file
#'
#' Post-hoc CRLF to LF conversion for files written by third-party functions
#' (e.g. \code{Biostrings::writeXStringSet()}) that may produce CRLF on Windows.
#' No-op on non-Windows platforms or if the file does not exist/is empty.
#'
#' @param path Character. File path to sanitise.
#'
#' @return Invisible \code{path}.
#' @keywords internal
ensure_lf <- function(path) {
  if (.Platform$OS.type == "windows" &&
      file.exists(path) &&
      file.size(path) > 0) {
    raw <- readBin(path, "raw", file.size(path))
    raw <- raw[raw != as.raw(0x0d)]
    writeBin(raw, path)
  }
  invisible(path)
}

#' Deduplicate protospacer sequences and build a guide index map
#'
#' Extracts unique protospacer sequences from a scored gRNA GRanges object
#' and builds a reverse index mapping each unique sequence to the GRanges
#' indices that share it. Used by off-target backends to search each unique
#' sequence once, then fan results back to all matching gRNAs.
#'
#' @param grna_gr GRanges. Must contain a \code{protospacer_sequence} metadata column.
#'
#' @return A named list:
#'   \describe{
#'     \item{all_spacers}{Character vector. All protospacers (uppercase),
#'           same length and order as \code{grna_gr}.}
#'     \item{unique_spacers}{Character vector. Deduplicated protospacers.}
#'     \item{guide_map}{Named list. Names = unique protospacers, values = integer
#'           vectors of GRanges indices sharing that sequence.}
#'   }
#' @keywords internal
deduplicate_protospacers <- function(grna_gr) {

  stopifnot(inherits(grna_gr, "GRanges"))
  stopifnot("protospacer_sequence" %in% names(GenomicRanges::mcols(grna_gr)))

  all_spacers <- toupper(as.character(GenomicRanges::mcols(grna_gr)$protospacer_sequence))
  unique_spacers <- unique(all_spacers)

  # Reverse index: unique sequence -> vector of GRanges positions
  guide_map <- lapply(
    setNames(unique_spacers, unique_spacers),
    function(seq) which(all_spacers == seq)
  )

  list(
    all_spacers    = all_spacers,
    unique_spacers = unique_spacers,
    guide_map      = guide_map
  )
}

#' Compute pair-level specificity using the additive off-target burden model
#'
#' For a two-guide deletion experiment, off-target damage from each gRNA is
#' independent and additive. The per-guide off-target burden is
#' \code{S = 100/specificity - 1} (the expected number of off-target cuts
#' relative to a perfect on-target cut). The pair burden is \code{S_pair = S_5p + S_3p},
#' and the pair specificity is \code{100 / (1 + S_pair)}.
#'
#' This is more stringent than \code{pmin()} for balanced pairs (correctly
#' reflecting doubled burden) and more lenient for asymmetric pairs (correctly
#' crediting a strong partner's low contribution to total risk).
#'
#' The formula implies an individual floor: a \code{min_pair_specificity} of T
#' requires each guide to have specificity >= \code{100 / (100/T - 1 + 1)} ~
#' \code{T} for small T. Guides with specificity near zero cannot form passing
#' pairs regardless of partner quality.
#'
#' @param spec_5p Numeric vector. 5-prime gRNA specificity scores (0-100 scale).
#' @param spec_3p Numeric vector. 3-prime gRNA specificity scores (0-100 scale).
#'
#' @return Numeric vector. Pair specificity on the same 0-100 MIT scale.
#'   Returns \code{NA} where either input is \code{NA} (off-target scoring
#'   not performed), and 0 where either input is <= 0.
#'
#' @examples
#' compute_pair_specificity(90, 90)   # 81.8 — two good guides
#' compute_pair_specificity(50, 50)   # 33.3 — doubled burden vs single
#' compute_pair_specificity(45, 90)   # 42.9 — strong partner compensates
#' compute_pair_specificity(30, 30)   # 17.6 — both mediocre
#'
#' @keywords internal
compute_pair_specificity <- function(spec_5p, spec_3p) {
  ifelse(is.na(spec_5p) | is.na(spec_3p),
         NA_real_,
         ifelse(spec_5p <= 0 | spec_3p <= 0,
                0,
                100 / (100 / spec_5p + 100 / spec_3p - 1)))
}

#' Convert a Windows file path to a WSL-compatible path
#'
#' Translates drive-letter paths (e.g. \code{C:/Users/...}) to WSL mount
#' paths (\code{/mnt/c/Users/...}). No-op if the path does not start with
#' a drive letter. Used by CRISPRitz backend functions that route through
#' WSL on Windows.
#'
#' @param win_path Character. A Windows-style file path.
#'
#' @return Character. The WSL-compatible equivalent.
#' @keywords internal
wsl_path <- function(win_path) {
  p <- normalizePath(win_path, winslash = "/", mustWork = FALSE)
  if (grepl("^[A-Za-z]:", p)) {
    drive <- tolower(substring(p, 1, 1))
    p <- paste0("/mnt/", drive, substring(p, 3))
  }
  p
}

#' Compute linguistic complexity of DNA sequences
#'
#' Calculates k-mer diversity as a proxy for sequence repetitiveness.
#' For each sequence, counts the number of unique k-mers (sliding window)
#' and divides by the total number of k-mer windows. Low values indicate
#' repetitive or low-complexity sequences (e.g. Alu-derived motifs).
#'
#' @param sequences Character vector. DNA sequences (uppercase recommended;
#'   mixed case is tolerated).
#' @param k Integer. K-mer size (default 3). Trigrams balance sensitivity
#'   (captures local repeats) with specificity (avoids flagging short
#'   motifs in otherwise complex sequences).
#'
#' @return Numeric vector, same length as \code{sequences}. Values range
#'   from \code{1/n_windows} (all k-mers identical, e.g. poly-T) to
#'   \code{1.0} (all k-mers unique). Returns \code{0} for sequences
#'   shorter than \code{k}.
#'
#' @examples
#' # Alu-derived repeat (low complexity)
#' compute_linguistic_complexity("TTTTTCTTTTTCTTTGAGAC")
#' # ~0.50 (9 unique trigrams / 18 windows)
#'
#' # Typical unique gRNA (high complexity)
#' compute_linguistic_complexity("AGCTTAGCGATCGTAGCTAG")
#' # ~0.94
#'
#' # Vectorised
#' compute_linguistic_complexity(c("TTTTTTTTTTTTTTTTTTTT", "AGCTTAGCGATCGTAGCTAG"))
#' # c(0.056, 0.94)
#'
#' @seealso \code{\link{screen_repeat_rich}} which uses this as one of two
#'   screening signals.
#'
#' @keywords internal
compute_linguistic_complexity <- function(sequences, k = 3L) {
  vapply(toupper(sequences), function(seq) {
    n <- nchar(seq)
    if (n < k) return(0)
    n_windows <- n - k + 1L
    starts <- seq_len(n_windows)
    kmers  <- substring(seq, starts, starts + k - 1L)
    length(unique(kmers)) / n_windows
  }, numeric(1), USE.NAMES = FALSE)
}


#' Screen protospacers for repeat-rich sequences likely to crash CRISPRitz
#'
#' Identifies guides that match highly repetitive genomic regions and would
#' cause CRISPRitz OOM during indexed bulge search (unbounded per-thread
#' result accumulation for guides matching millions of loci). Flagged guides
#' are quarantined with \code{specificity_score = 0} by the caller.
#'
#' Uses two complementary signals:
#' \describe{
#'   \item{Bowtie hit count}{(when available) Guides with more than
#'     \code{bowtie_threshold} off-target hits in mismatch-only search
#'     are in repetitive territory. Zero additional cost in hybrid mode
#'     since Phase 1 Bowtie data is already computed.}
#'   \item{Linguistic complexity}{Trigram diversity via
#'     \code{\link{compute_linguistic_complexity}}. Catches Alu-derived
#'     poly-T guides (LC ~0.50) without alignment data. Less sensitive
#'     to moderate-complexity repeat junctions (LC ~0.89); these require
#'     Bowtie counts or are caught reactively by bisection.}
#' }
#'
#' Either signal alone is sufficient to flag a guide.
#'
#' @param spacers Character vector. Protospacer sequences to screen.
#' @param bowtie_hit_counts Integer vector or NULL. Per-guide off-target
#'   hit counts from Bowtie Phase 1 (e.g. \code{n_offtargets} from
#'   \code{\link{score_offtargets}}). Same length as \code{spacers}.
#'   NULL when Bowtie data is unavailable (pure CRISPRitz mode).
#' @param bowtie_threshold Integer. Minimum Bowtie hit count to flag a
#'   guide as repeat-rich (default 5000). Typical on-target guides have
#'   1--100 hits at 3 mismatches; >5000 indicates repetitive sequence.
#' @param lc_threshold Numeric. Linguistic complexity below this value
#'   flags a guide (default 0.6). Calibrated to catch Alu-derived guides
#'   (LC 0.44--0.56) with margin while avoiding false positives on
#'   legitimate guides (LC > 0.85).
#' @param quiet Logical. Suppress diagnostic messages (default FALSE).
#'
#' @return Logical vector, same length as \code{spacers}. \code{TRUE}
#'   indicates a repeat-rich guide that should be quarantined.
#'
#' @seealso \code{\link{compute_linguistic_complexity}} for the complexity
#'   metric, \code{\link{bisect_failed_batch}} for reactive quarantine
#'   of guides that escape pre-screening.
#'
#' @keywords internal
screen_repeat_rich <- function(spacers,
                               bowtie_hit_counts = NULL,
                               bowtie_threshold = 5000L,
                               lc_threshold = 0.6,
                               quiet = FALSE) {

  n <- length(spacers)
  flagged_bt <- rep(FALSE, n)
  flagged_lc <- rep(FALSE, n)

  # ---- 1. Bowtie hit count signal (when available) ----
  if (!is.null(bowtie_hit_counts)) {
    if (length(bowtie_hit_counts) != n) {
      warning("bowtie_hit_counts length (", length(bowtie_hit_counts),
              ") does not match spacers length (", n, "). Ignoring Bowtie signal.")
    } else {
      flagged_bt <- !is.na(bowtie_hit_counts) & bowtie_hit_counts > bowtie_threshold
    }
  }

  # ---- 2. Linguistic complexity signal (always computed) ----
  lc <- compute_linguistic_complexity(spacers)
  flagged_lc <- lc < lc_threshold

  # ---- 3. Union of both signals ----
  flagged <- flagged_bt | flagged_lc

  # ---- 4. Diagnostic reporting ----
  if (!quiet && any(flagged)) {
    n_total <- sum(flagged)
    n_bt    <- sum(flagged_bt)
    n_lc    <- sum(flagged_lc)
    n_both  <- sum(flagged_bt & flagged_lc)

    detail_parts <- character(0)
    if (n_bt > 0)   detail_parts <- c(detail_parts,
                                      paste0(n_bt, " by Bowtie hit count >", bowtie_threshold))
    if (n_lc > 0)   detail_parts <- c(detail_parts,
                                      paste0(n_lc, " by linguistic complexity <", lc_threshold))
    if (n_both > 0)  detail_parts <- c(detail_parts,
                                       paste0(n_both, " by both signals"))

    message("Repeat-rich screen: ", n_total, "/", n, " guide(s) flagged (",
            paste(detail_parts, collapse = "; "), ")")
  }

  flagged
}
