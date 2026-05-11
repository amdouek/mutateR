#' @title Nuclease specification objects
#'
#' @description Construct, validate, and resolve `NucleaseSpec` objects that
#' describe a Cas effector for use with mutateR's site-finding and pair-assembly
#' steps. A spec captures the PAM/PFS topology, protospacer length, cut-site
#' geometry, and gRNA architecture sufficient to scan exonic sequence and
#' assemble candidate deletion pairs for any dsDNA cutter or nickase.
#'
#' Canonical specs (\code{"Cas9"}, \code{"Cas12a"}, \code{"enCas12a"}) preserve
#' existing pipeline behaviour, including on-target and off-target scoring.
#' User-supplied (custom) specs are intentionally restricted to pipeline stages
#' that do not require nuclease-specific trained models: site finding, phase
#' filtering, pair assembly, primer design, and plotting. On-target and
#' off-target scoring are skipped for custom specs.
#'
#' @name nuclease_spec
NULL


#' Construct a NucleaseSpec
#'
#' @param name Character. Display name (e.g. "AacCas12b", "SauCas9"). Required.
#' @param pam Character. PAM (or PFS) sequence with IUPAC ambiguity codes
#'        (e.g. "NGG", "TTTV", "NNGRRT", "TTN").
#' @param pam_side Character. Either \code{"3prime"} (PAM downstream of the
#'        protospacer on the target strand, Cas9-style) or \code{"5prime"}
#'        (PAM upstream, Cas12a-style).
#' @param protospacer_length Integer. Length of the protospacer in nucleotides.
#'        Must be in [10, 40].
#' @param cut_offset_top Integer. Signed offset of the top-strand cleavage
#'        position relative to the protospacer 3' end on the target strand
#'        (i.e., the protospacer-distal end). Negative values indicate cuts
#'        inside the protospacer (blunt Cas9 cleaves at -3); positive values
#'        indicate cuts downstream of the protospacer (Cas12a top strand cuts
#'        ~+18 from the PAM, equivalent to ~-5 from the protospacer 3' end
#'        for a 23-nt spacer).
#' @param cut_offset_bottom Integer. Signed offset of the bottom-strand
#'        cleavage position, in the same coordinate system as
#'        \code{cut_offset_top}. Equal to \code{cut_offset_top} for blunt
#'        cutters; different for staggered cutters (Cas12a, Cas12b).
#' @param target_type Character. One of \code{"dsDNA"}, \code{"ssDNA"},
#'        \code{"RNA"}. Only \code{"dsDNA"} is permitted for use with mutateR
#'        (exon-deletion pipeline). Defaults to \code{"dsDNA"}.
#' @param activity Character. One of \code{"cut"} (clean DSB), \code{"nick"}
#'        (single-strand nick — two nicks on opposite strands at flanking sites
#'        can substitute for a DSB), \code{"bind"} (catalytically dead), or
#'        \code{"unknown"}. Only \code{"cut"} and \code{"nick"} are permitted.
#' @param grna_architecture Character. Free-text description of the gRNA
#'        components needed (e.g. \code{"crRNA"}, \code{"crRNA + tracr"},
#'        \code{"wRNA"}). Stored for reporting only; does not affect logic.
#' @param pam_is_pfs Logical. \code{TRUE} if the sequence motif is a
#'        protospacer-flanking site (PFS) rather than a PAM (relevant for
#'        some RNA-targeting nucleases). Currently must be \code{FALSE};
#'        reserved for future use.
#' @param notes Character or NULL. Free-text annotation retained for reporting.
#' @param source Character or NULL. Provenance string (e.g. "CasPEDIA:8.4.2").
#'
#' @details
#' \strong{Cut-site convention.} Offsets are measured relative to the
#' \emph{protospacer 3' end} on the target strand. This unifies the two PAM
#' topologies into a single coordinate frame:
#' \itemize{
#'   \item \strong{3' PAM (Cas9):} target strand reads
#'         \code{[protospacer 5'] ... [protospacer 3'] [PAM]}.
#'         SpCas9 cuts the top strand 3 nt upstream of the PAM (i.e. inside the
#'         protospacer) and the bottom strand at the same position
#'         (\code{cut_offset_top = cut_offset_bottom = -3}). Blunt cut.
#'   \item \strong{5' PAM (Cas12a):} target strand reads
#'         \code{[PAM] [protospacer 5'] ... [protospacer 3']}.
#'         AsCas12a cuts the top strand ~5 nt before the protospacer 3' end
#'         (\code{cut_offset_top = -5}) and the bottom strand ~9 nt before
#'         (\code{cut_offset_bottom = -9}). Staggered cut leaving a 4-nt 5'
#'         overhang.
#' }
#' The single \code{cut_site} integer used by downstream code is the midpoint
#' \code{floor((top + bottom) / 2)}. The finder also emits genomic
#' \code{cut_site_top} and \code{cut_site_bottom} columns so the full stagger
#' geometry remains recoverable.
#'
#' \strong{Allowed activity types.} Custom specs are restricted to
#' \code{"cut"} and \code{"nick"} because mutateR's downstream logic
#' (deletion-pair assembly, primer design) is meaningful only for
#' nucleases that produce double-strand breaks (or paired nicks). RNA-targeting
#' and binding-only effectors are rejected at construction.
#'
#' @return An object of class \code{"NucleaseSpec"} (a structured list).
#'
#' @examples
#' \dontrun{
#' # SauCas9 (CasPEDIA 1.1.1): 3' NNGRRT PAM, 21-nt spacer, blunt cut at -3
#' sau <- nuclease_spec(
#'   name               = "SauCas9",
#'   pam                = "NNGRRT",
#'   pam_side           = "3prime",
#'   protospacer_length = 21L,
#'   cut_offset_top     = -3L,
#'   cut_offset_bottom  = -3L,
#'   grna_architecture  = "crRNA + tracr",
#'   source             = "CasPEDIA:1.1.1"
#' )
#'
#' # AacCas12b (CasPEDIA 8.4.2): 5' TTN PAM, 20-nt spacer, staggered cut
#' aac <- nuclease_spec(
#'   name               = "AacCas12b",
#'   pam                = "TTN",
#'   pam_side           = "5prime",
#'   protospacer_length = 20L,
#'   cut_offset_top     = -5L,
#'   cut_offset_bottom  = -9L,
#'   activity           = "cut",
#'   grna_architecture  = "crRNA + tracr",
#'   source             = "CasPEDIA:8.4.2"
#' )
#' }
#'
#' @export
nuclease_spec <- function(name,
                          pam,
                          pam_side,
                          protospacer_length,
                          cut_offset_top,
                          cut_offset_bottom,
                          target_type       = "dsDNA",
                          activity          = "cut",
                          grna_architecture = "crRNA",
                          pam_is_pfs        = FALSE,
                          notes             = NULL,
                          source            = NULL) {

  # ---- Required-field presence ----
  if (missing(name) || !is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.")
  }
  if (missing(pam) || !is.character(pam) || length(pam) != 1L || !nzchar(pam)) {
    stop("`pam` must be a single non-empty character string (IUPAC DNA).")
  }
  if (missing(pam_side) || !pam_side %in% c("3prime", "5prime")) {
    stop("`pam_side` must be either '3prime' or '5prime'.")
  }
  if (missing(protospacer_length) || !is_single_integerish(protospacer_length)) {
    stop("`protospacer_length` must be a single integer.")
  }
  if (missing(cut_offset_top) || !is_single_integerish(cut_offset_top)) {
    stop("`cut_offset_top` must be a single integer.")
  }
  if (missing(cut_offset_bottom) || !is_single_integerish(cut_offset_bottom)) {
    stop("`cut_offset_bottom` must be a single integer.")
  }

  # ---- Value-domain checks ----
  pam_upper <- toupper(pam)
  if (!grepl("^[ACGTRYSWKMBDHVN]+$", pam_upper)) {
    stop("`pam` must contain only IUPAC DNA characters (ACGT plus ambiguity codes).")
  }

  protospacer_length <- as.integer(protospacer_length)
  if (protospacer_length < 10L || protospacer_length > 40L) {
    stop("`protospacer_length` must be between 10 and 40 (got ", protospacer_length, ").")
  }

  cut_offset_top    <- as.integer(cut_offset_top)
  cut_offset_bottom <- as.integer(cut_offset_bottom)

  # ---- Capability gating (target_type and activity) ----
  if (!target_type %in% c("dsDNA", "ssDNA", "RNA")) {
    stop("`target_type` must be one of 'dsDNA', 'ssDNA', 'RNA' (got '", target_type, "').")
  }
  if (target_type != "dsDNA") {
    stop("mutateR currently supports dsDNA-targeting nucleases only. ",
         "Got target_type = '", target_type, "'. ",
         "RNA- and ssDNA-targeting effectors are out of scope for ",
         "the exon-deletion pipeline.")
  }

  if (!activity %in% c("cut", "nick", "bind", "unknown")) {
    stop("`activity` must be one of 'cut', 'nick', 'bind', 'unknown' (got '", activity, "').")
  }
  if (!activity %in% c("cut", "nick")) {
    stop("Custom NucleaseSpec requires `activity` to be 'cut' or 'nick'. ",
         "Got '", activity, "'. Binding-only and unknown-activity effectors ",
         "cannot drive deletion-pair assembly and are rejected at construction.")
  }

  if (!isTRUE(pam_is_pfs) && !identical(pam_is_pfs, FALSE)) {
    stop("`pam_is_pfs` must be TRUE or FALSE.")
  }
  if (isTRUE(pam_is_pfs)) {
    stop("`pam_is_pfs = TRUE` is not supported in this release ",
         "(reserved for PFS-style RNA-targeting nucleases).")
  }

  if (!is.character(grna_architecture) || length(grna_architecture) != 1L) {
    stop("`grna_architecture` must be a single character string.")
  }
  if (!is.null(notes)  && (!is.character(notes)  || length(notes)  != 1L)) {
    stop("`notes` must be NULL or a single character string.")
  }
  if (!is.null(source) && (!is.character(source) || length(source) != 1L)) {
    stop("`source` must be NULL or a single character string.")
  }

  spec <- list(
    name               = name,
    pam                = pam_upper,
    pam_side           = pam_side,
    pam_length         = nchar(pam_upper),
    protospacer_length = protospacer_length,
    cut_offset_top     = cut_offset_top,
    cut_offset_bottom  = cut_offset_bottom,
    target_type        = target_type,
    activity           = activity,
    grna_architecture  = grna_architecture,
    pam_is_pfs         = FALSE,
    notes              = notes,
    source             = source,
    is_canonical       = FALSE,
    canonical_key      = NA_character_
  )
  class(spec) <- c("NucleaseSpec", "list")
  spec
}


#' @keywords internal
is_single_integerish <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && x == as.integer(x)
}


#' Built-in NucleaseSpec for SpCas9 (canonical mutateR "Cas9" path)
#' @keywords internal
.spec_Cas9 <- function() {
  s <- nuclease_spec(
    name               = "SpCas9",
    pam                = "NGG",
    pam_side           = "3prime",
    protospacer_length = 20L,
    cut_offset_top     = -3L,
    cut_offset_bottom  = -3L,
    activity           = "cut",
    grna_architecture  = "crRNA + tracr",
    source             = "canonical:Cas9"
  )
  s$is_canonical  <- TRUE
  s$canonical_key <- "Cas9"
  s
}

#' Built-in NucleaseSpec for AsCas12a (canonical mutateR "Cas12a" path)
#'
#' @details AsCas12a top strand cuts ~18 nt downstream of the PAM (between
#' positions 18 and 19 from PAM end); bottom strand cuts ~23 nt downstream
#' (between positions 22 and 23). For a 23-nt protospacer (PAM end = +0,
#' protospacer end = +23), this translates to offsets relative to the
#' protospacer 3' end of approximately -5 (top) and -1 (bottom).
#' @keywords internal
.spec_Cas12a <- function() {
  s <- nuclease_spec(
    name               = "AsCas12a",
    pam                = "TTTV",
    pam_side           = "5prime",
    protospacer_length = 23L,
    cut_offset_top     = -5L,
    cut_offset_bottom  = -1L,
    activity           = "cut",
    grna_architecture  = "crRNA",
    source             = "canonical:Cas12a"
  )
  s$is_canonical  <- TRUE
  s$canonical_key <- "Cas12a"
  s
}

#' Built-in NucleaseSpec for enAsCas12a (canonical mutateR "enCas12a" path)
#' @keywords internal
.spec_enCas12a <- function() {
  s <- nuclease_spec(
    name               = "enAsCas12a",
    pam                = "TTTN",
    pam_side           = "5prime",
    protospacer_length = 23L,
    cut_offset_top     = -5L,
    cut_offset_bottom  = -1L,
    activity           = "cut",
    grna_architecture  = "crRNA",
    source             = "canonical:enCas12a"
  )
  s$is_canonical  <- TRUE
  s$canonical_key <- "enCas12a"
  s
}


#' Resolve a user-facing `nuclease` argument to a NucleaseSpec
#'
#' Accepts the canonical string keys (\code{"Cas9"}, \code{"Cas12a"},
#' \code{"enCas12a"}) or a pre-built \code{NucleaseSpec} object and returns
#' the corresponding spec. Used internally by site finders and pipeline
#' entry points to normalise the nuclease argument.
#'
#' @param nuclease Character string (one of "Cas9", "Cas12a", "enCas12a") or
#'        a \code{NucleaseSpec} object.
#'
#' @return A \code{NucleaseSpec} object.
#' @keywords internal
resolve_nuclease <- function(nuclease) {
  if (inherits(nuclease, "NucleaseSpec")) {
    return(nuclease)
  }
  if (is.character(nuclease) && length(nuclease) == 1L) {
    return(switch(nuclease,
                  "Cas9"     = .spec_Cas9(),
                  "Cas12a"   = .spec_Cas12a(),
                  "enCas12a" = .spec_enCas12a(),
                  stop("Unknown nuclease string '", nuclease, "'. ",
                       "Expected one of 'Cas9', 'Cas12a', 'enCas12a', ",
                       "or a NucleaseSpec object from nuclease_spec().")))
  }
  stop("`nuclease` must be a character string ('Cas9', 'Cas12a', 'enCas12a') ",
       "or a NucleaseSpec object.")
}


#' Print method for NucleaseSpec
#' @param x A NucleaseSpec object.
#' @param ... Unused.
#' @export
print.NucleaseSpec <- function(x, ...) {
  cat("<NucleaseSpec>\n")
  cat("  Name              : ", x$name, "\n", sep = "")
  cat("  PAM               : ", x$pam, " (", x$pam_side, ")\n", sep = "")
  cat("  Protospacer length: ", x$protospacer_length, " nt\n", sep = "")
  if (x$cut_offset_top == x$cut_offset_bottom) {
    cat("  Cleavage          : blunt, offset ",
        x$cut_offset_top, " from protospacer 3' end\n", sep = "")
  } else {
    cat("  Cleavage          : staggered (top ", x$cut_offset_top,
        ", bottom ", x$cut_offset_bottom,
        " from protospacer 3' end)\n", sep = "")
  }
  cat("  Activity          : ", x$activity, "\n", sep = "")
  cat("  Target            : ", x$target_type, "\n", sep = "")
  cat("  gRNA architecture : ", x$grna_architecture, "\n", sep = "")
  if (!is.null(x$source)) cat("  Source            : ", x$source, "\n", sep = "")
  if (!is.null(x$notes))  cat("  Notes             : ", x$notes,  "\n", sep = "")
  cat("  Canonical         : ", if (isTRUE(x$is_canonical)) "yes" else "no", "\n", sep = "")
  invisible(x)
}
