#' @title Genome cache management for off-target scoring
#'
#' @description Prepares and persistently caches genome data for off-target
#' search backends. Currently provides BSgenome to per-chromosome FASTA export
#' (consumed by CRISPRitz and used as source for Bowtie index building).
#'
#' @keywords internal
#' @name ot_cache
NULL


#' Export BSgenome to per-chromosome FASTA files (cached)
#'
#' Exports standard chromosomes from a BSgenome object to individual FASTA files
#' in a persistent cache directory. Subsequent calls with the same genome return
#' the cached path immediately. A sentinel file records cache completeness.
#'
#' @param genome BSgenome object.
#' @param cache_dir Character. Base cache directory.
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return Character. Path to the directory containing per-chromosome FASTA files.
#' @keywords internal
ensure_genome_cache <- function(genome, cache_dir, quiet = FALSE) {

  stopifnot(inherits(genome, "BSgenome"))

  # ---- 1. Derive cache subdirectory from BSgenome identity ----
  genome_name <- genome@pkgname  # e.g. "BSgenome.Hsapiens.UCSC.hg38"
  genome_dir  <- file.path(cache_dir, genome_name)
  sentinel    <- file.path(cache_dir, paste0(genome_name, ".complete"))

  # ---- 2. Check for valid existing cache ----
  if (file.exists(sentinel)) {
    meta <- tryCatch(
      readLines(sentinel, warn = FALSE),
      error = function(e) character(0)
    )

    if (length(meta) >= 1 && meta[1] == genome_name) {
      if (!quiet) message("Using cached genome FASTA: ", genome_dir)
      return(genome_dir)
    } else {
      if (!quiet) message("Cache metadata mismatch. Rebuilding genome FASTA cache...")
      unlink(genome_dir, recursive = TRUE)
    }
  }

  # ---- 3. Create output directory ----
  if (!dir.exists(genome_dir)) {
    dir.create(genome_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- 4. Identify standard chromosomes ----
  std_chroms <- GenomeInfoDb::standardChromosomes(genome)

  if (length(std_chroms) == 0) {
    stop("No standard chromosomes found in genome object '", genome_name,
         "'. Check that the BSgenome object is valid.")
  }

  if (!quiet) message("Exporting ", length(std_chroms),
                      " standard chromosomes from ", genome_name, " to FASTA...")

  # ---- 5. Export each chromosome as individual FASTA ----
  for (chrom in std_chroms) {

    fasta_path <- file.path(genome_dir, paste0(chrom, ".fa"))

    # Skip if this chromosome was already exported (partial cache recovery)
    if (file.exists(fasta_path) && file.size(fasta_path) > 0) {
      next
    }

    chrom_seq <- tryCatch(
      Biostrings::getSeq(genome, chrom),
      error = function(e) {
        warning("Failed to extract sequence for ", chrom, ": ", e$message)
        return(NULL)
      }
    )

    if (is.null(chrom_seq)) next

    chrom_set <- Biostrings::DNAStringSet(chrom_seq)
    names(chrom_set) <- chrom

    Biostrings::writeXStringSet(chrom_set, filepath = fasta_path, format = "fasta")
    ensure_lf(fasta_path)
  }

  # ---- 6. Write chromosome lengths file ----
  # Stored OUTSIDE the FASTA directory — CRISPRitz rejects non-.fa files
  chrom_lengths <- GenomeInfoDb::seqlengths(genome)[std_chroms]
  lengths_file  <- file.path(cache_dir, paste0(genome_name, ".chrom.sizes"))
  lengths_lines <- paste(names(chrom_lengths), as.integer(chrom_lengths), sep = "\t")
  write_unix_lines(lengths_lines, lengths_file)

  # ---- 7. Write sentinel file ----
  # Stored OUTSIDE the FASTA directory — CRISPRitz rejects non-.fa files
  write_unix_lines(c(genome_name, as.character(Sys.time())), sentinel)

  if (!quiet) message("Genome FASTA cache complete: ", genome_dir)

  return(genome_dir)
}

#' Check whether a complete Bowtie index exists at a given prefix
#'
#' Bowtie creates six index files per genome. For genomes exceeding 4GB,
#' the large-index format (\code{.ebwtl}) is used instead of \code{.ebwt}.
#' This helper checks for a complete set in either format.
#'
#' @param prefix Character. Bowtie index prefix (path without suffix).
#'
#' @return Logical. TRUE if a complete set of index files exists.
#' @keywords internal
bowtie_index_exists <- function(prefix) {
  small_suffixes <- c(".1.ebwt",  ".2.ebwt",  ".3.ebwt",  ".4.ebwt",
                      ".rev.1.ebwt",  ".rev.2.ebwt")
  large_suffixes <- c(".1.ebwtl", ".2.ebwtl", ".3.ebwtl", ".4.ebwtl",
                      ".rev.1.ebwtl", ".rev.2.ebwtl")

  all(file.exists(paste0(prefix, small_suffixes))) ||
    all(file.exists(paste0(prefix, large_suffixes)))
}


#' Build and cache a Bowtie genome index for off-target searching
#'
#' Creates a Bowtie index from the per-chromosome FASTA files produced by
#' \code{\link{ensure_genome_cache}}. The index is built once and cached
#' persistently. Subsequent calls with the same genome return the cached
#' index prefix immediately.
#'
#' The returned prefix is suitable for passing directly to
#' \code{crisprBowtie::runCrisprBowtie(bowtie_index = ...)}.
#'
#' @param genome BSgenome object.
#' @param cache_dir Character. Base cache directory (same as used by
#'        \code{ensure_genome_cache}).
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return Character. Path to the Bowtie index prefix (i.e. the common
#'   prefix of the six \code{.ebwt}/\code{.ebwtl} files). Pass this
#'   directly as the \code{bowtie_index} argument to
#'   \code{crisprBowtie::runCrisprBowtie()}.
#'
#' @details
#' \strong{Build time:} Expected 10–15 minutes for the human genome
#' (hg38, 25 standard chromosomes) on Linux/macOS. Up to 2-3 hours on Windows.
#' Runs once per genome; the index persists across R sessions.
#'
#' \strong{Disk usage:} approximately 2.5 GB for hg38 (small-index format).
#' Genomes exceeding 4 GB total sequence (e.g. axolotl) use Bowtie's
#' large-index format (\code{.ebwtl}) automatically.
#'
#' \strong{Dependencies:} Requires the \pkg{Rbowtie} package
#' (Bioconductor). Checked at runtime with an informative error message.
#'
#' @seealso \code{\link{ensure_genome_cache}} for FASTA export,
#'          \code{\link[Rbowtie]{bowtie_build}} for the underlying builder.
#'
#' @keywords internal
ensure_bowtie_index <- function(genome, cache_dir, quiet = FALSE) {

  stopifnot(inherits(genome, "BSgenome"))

  # ---- 1. Check Rbowtie availability ----
  if (!requireNamespace("Rbowtie", quietly = TRUE)) {
    stop("Package 'Rbowtie' is required for the Bowtie off-target backend.\n",
         "Install with: BiocManager::install('Rbowtie')")
  }

  # ---- 2. Ensure FASTA cache exists ----
  # ensure_genome_cache() is idempotent — returns immediately if cache is valid
  genome_dir <- ensure_genome_cache(genome, cache_dir, quiet = quiet)

  # ---- 3. Derive index paths ----
  # Note: the actual index prefix depends on Rbowtie's naming convention
  # (derived from the first reference filename). It is read from the sentinel
  # on cache hit, or detected from build output on cache miss.
  genome_name <- genome@pkgname
  index_dir   <- file.path(cache_dir, paste0(genome_name, "_bowtie_index"))
  sentinel    <- file.path(index_dir, ".bowtie_index.complete")

  # ---- 4. Check for valid existing index ----
  if (file.exists(sentinel)) {
    meta <- tryCatch(
      readLines(sentinel, warn = FALSE),
      error = function(e) character(0)
    )

    # Sentinel format: line 1 = genome_name, line 2 = index_prefix, line 3 = timestamp
    if (length(meta) >= 2 && meta[1] == genome_name) {
      stored_prefix <- meta[2]
      if (bowtie_index_exists(stored_prefix)) {
        if (!quiet) message("Using cached Bowtie index: ", stored_prefix)
        return(stored_prefix)
      } else {
        if (!quiet) message("Bowtie index sentinel found but index files are missing. Rebuilding...")
        unlink(sentinel)
      }
    } else {
      if (!quiet) message("Bowtie index metadata mismatch. Rebuilding...")
      unlink(index_dir, recursive = TRUE)
    }
  }

  # ---- 5. Create index directory ----
  if (!dir.exists(index_dir)) {
    dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- 6. Collect FASTA files from genome cache ----
  fasta_files <- sort(list.files(genome_dir, pattern = "\\.fa$", full.names = TRUE))

  if (length(fasta_files) == 0) {
    stop("No FASTA files found in genome cache: ", genome_dir, "\n",
         "The FASTA cache may be corrupted. Try deleting\n  ",
         genome_dir, "\nand re-running to rebuild it.")
  }

  # ---- 7. Build index ----
  if (!quiet) {
    message("Building Bowtie index from ", length(fasta_files),
            " chromosomes in ", genome_name, "...")
    message("This is a one-time operation (~10-15 min for human genome). ",
            "Index will be cached for future sessions.")
  }

  t_start <- Sys.time()

  build_status <- tryCatch(
    Rbowtie::bowtie_build(
      references = fasta_files,
      outdir     = index_dir,
      force      = TRUE
    ),
    error = function(e) {
      stop("Bowtie index build failed: ", e$message, "\n",
           "Ensure the Rbowtie package is properly installed and ",
           "sufficient disk space is available (~2.5 GB for human genome).")
    }
  )

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  if (!quiet) message("Bowtie index built in ", round(as.numeric(t_elapsed), 1), " minutes.")

  # ---- 8. Detect and validate actual index prefix ----
  # Rbowtie derives the index prefix from the first reference filename:
  #   file.path(outdir, sub("\\.[^.]*$", "", basename(references[1])))
  # Rather than replicating that logic, detect the prefix from the generated files.
  ebwt_hits <- list.files(index_dir, pattern = "\\.1\\.ebwt[l]?$", full.names = TRUE)

  if (length(ebwt_hits) == 0) {
    created_files <- list.files(index_dir, full.names = FALSE)
    stop("Bowtie index build completed but no index files found.\n",
         "Index directory: ", index_dir, "\n",
         "Directory contents: ",
         if (length(created_files) > 0) paste(created_files, collapse = ", ") else "(empty)")
  }

  # Strip the .1.ebwt or .1.ebwtl suffix to recover the prefix
  index_prefix <- sub("\\.1\\.ebwt[l]?$", "", ebwt_hits[1])

  if (!bowtie_index_exists(index_prefix)) {
    created_files <- list.files(index_dir, full.names = FALSE)
    stop("Bowtie index build produced partial output.\n",
         "Detected prefix: ", index_prefix, "\n",
         "Directory contents: ", paste(created_files, collapse = ", "))
  }

  # ---- 9. Write sentinel ----
  # Sentinel format: line 1 = genome_name, line 2 = index_prefix, line 3 = timestamp
  write_unix_lines(c(genome_name, index_prefix, as.character(Sys.time())), sentinel)

  if (!quiet) message("Bowtie index cached at: ", index_prefix)

  return(index_prefix)
}

#' Build and cache a CRISPRitz genome index for bulge-aware off-target searching
#'
#' Pre-indexes a genome for CRISPRitz bulge-mode search, dramatically
#' accelerating subsequent off-target searches at the cost of a one-time
#' indexing step and additional disk usage. The index is PAM-specific
#' (one per nuclease type) and bulge-size-specific. Results are cached
#' persistently with sentinel-based validation.
#'
#' When \code{max_bulge = 0} (mismatch-only mode), indexing is skipped
#' entirely — CRISPRitz handles mismatch-only search efficiently without
#' pre-indexing.
#'
#' @param genome BSgenome object.
#' @param nuclease Character. One of \code{"Cas9"}, \code{"Cas12a"},
#'   \code{"enCas12a"}. Determines the PAM specification for indexing.
#' @param max_bulge Integer. Maximum bulge size to index for (default 2).
#'   Must match or exceed the \code{max_dna_bulges} / \code{max_rna_bulges}
#'   used in the subsequent search. Set to 0 to skip indexing.
#' @param cache_dir Character. Base cache directory (same as used by
#'   \code{\link{ensure_genome_cache}}).
#' @param threads Integer. Thread count for CRISPRitz indexing (default 4).
#' @param timeout Integer. Maximum seconds for indexing (default 3600).
#'   Genome indexing is slower than searching; ~20-60 min for hg38 is typical.
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return A named list:
#'   \describe{
#'     \item{index_dir}{Character. Path to the directory to pass to
#'           \code{run_crispritz_search()}. If indexing succeeded, this is
#'           the enriched/indexed directory. If indexing was skipped or
#'           failed, this is the raw FASTA directory (graceful degradation).}
#'     \item{indexed}{Logical. Whether a valid pre-built index is available.}
#'   }
#'
#' @details
#' \strong{Disk usage.} The CRISPRitz enriched genome index is typically
#' 2-4 GB for hg38 (depending on bulge size), stored alongside the existing
#' FASTA cache (~3 GB) and Bowtie index (~2.9 GB). Total cache footprint
#' with all backends: ~9-10 GB for human genome.
#'
#' \strong{PAM specificity.} Unlike the Bowtie index (which is genome-universal),
#' CRISPRitz indexes are PAM-specific. Switching nuclease types requires a
#' separate index. Each is cached independently with its own sentinel.
#'
#' @seealso \code{\link{ensure_genome_cache}} for FASTA export,
#'          \code{\link{ensure_bowtie_index}} for the analogous Bowtie function,
#'          \code{\link{run_crispritz_search}} for the search step.
#'
#' @keywords internal
ensure_crispritz_index <- function(genome,
                                   nuclease = c("Cas9", "Cas12a", "enCas12a"),
                                   max_bulge = 2L,
                                   cache_dir,
                                   threads = 4L,
                                   timeout = 3600L,
                                   quiet = FALSE) {

  stopifnot(inherits(genome, "BSgenome"))
  nuclease <- match.arg(nuclease)

  # ---- 1. Ensure FASTA cache ----
  genome_dir  <- ensure_genome_cache(genome, cache_dir, quiet = quiet)
  genome_name <- genome@pkgname

  # ---- 2. Skip if mismatch-only ----
  if (max_bulge <= 0L) {
    if (!quiet) message("max_bulge = 0: CRISPRitz genome indexing not required.")
    return(list(index_dir = genome_dir, indexed = FALSE))
  }

  # ---- 3. Check for valid cached index via sentinel ----
  index_id <- paste0(genome_name, "_crispritz_", nuclease, "_b", max_bulge)
  sentinel <- file.path(cache_dir, paste0(index_id, ".complete"))

  if (file.exists(sentinel)) {
    meta <- tryCatch(
      readLines(sentinel, warn = FALSE),
      error = function(e) character(0)
    )

    # Sentinel format: line 1 = index_id, line 2 = index_dir (TST subdir), line 3 = timestamp
    if (length(meta) >= 2 && meta[1] == index_id) {
      stored_dir <- meta[2]

      if (dir.exists(stored_dir)) {
        # ---- Reverse migration: fix sentinels erroneously pointing to genome_library/ parent ----
        # CRISPRitz parses the directory basename to extract PAM and max_bulges
        # (convention: {PAM}_{maxBulges}_{name}), so the sentinel must store the
        # TST subdirectory, not the genome_library/ parent.
        if (basename(stored_dir) == "genome_library") {
          pam_label <- switch(nuclease,
                              "Cas9"     = "NGG",
                              "Cas12a"   = "TTTV",
                              "enCas12a" = "TTTN")
          tst_name <- paste0(pam_label, "_", max_bulge, "_", index_id)
          tst_candidate <- file.path(stored_dir, tst_name)

          if (dir.exists(tst_candidate) &&
              length(list.files(tst_candidate, pattern = "\\.bin$")) > 0) {
            if (!quiet) message("Correcting CRISPRitz index sentinel: ",
                                "genome_library/ -> ", tst_name)
            write_unix_lines(c(index_id, tst_candidate, as.character(Sys.time())), sentinel)
            return(list(index_dir = tst_candidate, indexed = TRUE))
          } else {
            if (!quiet) message("Sentinel points to genome_library/ but expected TST subdir\n  ",
                                tst_candidate, "\n  not found. Will attempt rebuild.")
            unlink(sentinel)
          }
        } else {
          # Normal case: sentinel points to TST subdir — use as-is
          if (!quiet) message("Using cached CRISPRitz index: ", stored_dir)
          return(list(index_dir = stored_dir, indexed = TRUE))
        }

      } else {
        if (!quiet) message("CRISPRitz index sentinel points to non-existent directory:\n  ",
                            stored_dir, "\nWill attempt to rebuild index.")
        unlink(sentinel)
      }

    } else {
      if (!quiet) {
        if (length(meta) < 2) {
          message("CRISPRitz index sentinel is malformed (", length(meta),
                  " lines). Rebuilding...")
        } else {
          message("CRISPRitz index sentinel ID mismatch (expected '", index_id,
                  "', found '", meta[1], "'). Rebuilding...")
        }
      }
      unlink(sentinel)
    }
  }

  # ---- 4. Locate CRISPRitz executable ----
  cz <- locate_crispritz_bin(quiet = quiet)

  if (is.null(cz$bin)) {
    warning("CRISPRitz executable not found. Cannot pre-index genome.\n",
            "Bulge-mode search will proceed without pre-indexing (significantly slower).\n",
            "To enable indexing, install CRISPRitz via install_mutater_env().")
    return(list(index_dir = genome_dir, indexed = FALSE))
  }

  if (!quiet) message("Using CRISPRitz: ", cz$bin,
                      if (cz$use_wsl) " (via WSL)" else "")

  # ---- 5. Write PAM file for this nuclease ----
  guide_len <- switch(nuclease,
                      "Cas9"     = 20L,
                      "Cas12a"   = 23L,
                      "enCas12a" = 23L)

  guide_placeholder <- paste(rep("N", guide_len), collapse = "")

  pam_len <- switch(nuclease,
                    "Cas9"     = 3L,
                    "Cas12a"   = 4L,
                    "enCas12a" = 4L)

  pam_string <- switch(nuclease,
                       "Cas9"     = paste0(guide_placeholder, "NGG"),
                       "Cas12a"   = paste0("TTTV", guide_placeholder),
                       "enCas12a" = paste0("TTTN", guide_placeholder))

  # CRISPRitz prefixes output with {PAM}_{bMax}_{name} — capture PAM label
  pam_label <- switch(nuclease,
                      "Cas9"     = "NGG",
                      "Cas12a"   = "TTTV",
                      "enCas12a" = "TTTN")

  tmpdir <- file.path(tempdir(),
                      paste0("mutateR_cz_idx_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

  pam_file <- file.path(tmpdir, "pam.txt")
  write_unix_lines(paste(pam_string, pam_len), pam_file)

  # ---- 6. Build index command ----
  if (!quiet) {
    message("Building CRISPRitz genome index for ", nuclease,
            " (bMax=", max_bulge, ")...")
    message("This is a one-time operation per nuclease type. ",
            "Expected time: ~20-60 min for human genome.")
    message("Expected disk usage: ~2-4 GB additional.")
  }

  # CRISPRitz index-genome requires: name_genome genomeDir pamFile -bMax N
  # name_genome is a bare label (no path separators!) used to name output.
  # CRISPRitz creates output in CWD as:
  #   {name_genome}_enriched/           — enriched FASTA files
  #   {PAM}_{bMax}_{name_genome}/       — TST index files
  # We cd to cache_dir before execution so all output lands there.
  index_name <- index_id

  args <- c("index-genome", index_name, genome_dir, pam_file,
            "-bMax", as.character(max_bulge))

  stdout_log <- file.path(tmpdir, "index_stdout.log")
  stderr_log <- file.path(tmpdir, "index_stderr.log")

  t_start <- Sys.time()

  # ---- 7. Execute ----
  if (cz$use_wsl) {

    wsl_cache_dir <- wsl_path(cache_dir)

    wsl_args <- c("index-genome",
                  index_name,
                  wsl_path(genome_dir),
                  wsl_path(pam_file),
                  "-bMax", as.character(max_bulge))

    wsl_cmd_str <- paste(
      "source ~/miniconda3/etc/profile.d/conda.sh",
      "&&", "conda activate r-mutater-wsl",
      "&&", "cd", wsl_cache_dir,
      "&&", "timeout", as.character(timeout),
      "crispritz.py", paste(wsl_args, collapse = " "),
      ">", wsl_path(stdout_log),
      "2>", wsl_path(stderr_log)
    )

    result <- system2("wsl",
                      args = c("bash", "-lc", shQuote(wsl_cmd_str)),
                      stdout = TRUE, stderr = TRUE)
    exit_code <- attr(result, "status")
    if (is.null(exit_code)) exit_code <- 0L

  } else {

    # cd to cache_dir so CRISPRitz creates output there
    old_wd <- setwd(cache_dir)
    on.exit(setwd(old_wd), add = TRUE)

    if (grepl("\\.py$", cz$bin)) {
      exec      <- cz$python
      full_args <- c(cz$bin, args)
    } else {
      exec      <- cz$bin
      full_args <- args
    }

    exit_code <- system2(exec, full_args,
                         stdout = stdout_log, stderr = stderr_log)
  }

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")

  # ---- 8. Diagnostic logging ----
  if (!quiet) {
    for (logfile in c(stdout_log, stderr_log)) {
      if (file.exists(logfile)) {
        content <- tryCatch(readLines(logfile, warn = FALSE, n = 10),
                            error = function(e) character(0))
        if (length(content) > 0 && any(nzchar(content))) {
          label <- if (grepl("stdout", logfile)) "stdout" else "stderr"
          message("CRISPRitz index-genome ", label, ":\n  ",
                  paste(content, collapse = "\n  "))
        }
      }
    }
  }

  # Temp diagnostic: Log cache_dir contents immediately after execution (before any cleanup)
  if (!quiet) {
    post_contents <- list.dirs(cache_dir, recursive = FALSE, full.names = FALSE)
    message("Cache directory contents after indexing: ", paste(post_contents, collapse = ", "))
  }

  # ---- 9. Check execution status ----
  if (exit_code == 124) {
    warning("CRISPRitz genome indexing timed out after ", timeout, " seconds.\n",
            "Consider increasing timeout or indexing on a faster machine.\n",
            "Search will proceed without pre-indexing (slower).")
    return(list(index_dir = genome_dir, indexed = FALSE))
  }

  if (exit_code != 0) {
    err_msg <- tryCatch(
      paste(readLines(stderr_log, warn = FALSE), collapse = "\n"),
      error = function(e) "(could not read log)")
    warning("CRISPRitz genome indexing failed (exit code ", exit_code, "):\n",
            err_msg, "\nSearch will proceed without pre-indexing (slower).")
    return(list(index_dir = genome_dir, indexed = FALSE))
  }

  if (!quiet) message("CRISPRitz index-genome completed in ",
                      round(as.numeric(t_elapsed), 1), " minutes.")

  # ---- 10. Detect output directory ----
  # CRISPRitz search -index expects the TST index directory directly.
  # It parses basename to extract PAM and max_bulges (convention:
  # {PAM}_{bMax}_{name}). We must store and return the TST subdir path.
  genome_lib <- file.path(cache_dir, "genome_library")
  tst_name   <- paste0(pam_label, "_", max_bulge, "_", index_name)
  tst_dir    <- file.path(genome_lib, tst_name)

  index_dir <- NULL

  # Primary check: expected TST subdirectory with .bin files
  if (dir.exists(tst_dir)) {
    bin_files <- list.files(tst_dir, pattern = "\\.bin$", recursive = FALSE)
    if (length(bin_files) > 0) {
      index_dir <- tst_dir
      if (!quiet) message("Found TST index: ", tst_dir,
                          " (", length(bin_files), " .bin files)")
    }
  }

  # Fallback: scan genome_library/ subdirectories for any with .bin files
  if (is.null(index_dir) && dir.exists(genome_lib)) {
    gl_subdirs <- list.dirs(genome_lib, recursive = FALSE, full.names = TRUE)
    for (subdir in gl_subdirs) {
      bin_files <- list.files(subdir, pattern = "\\.bin$", recursive = FALSE)
      if (length(bin_files) > 0) {
        index_dir <- subdir
        if (!quiet) message("Found index files in: ", basename(subdir),
                            " (", length(bin_files), " .bin files)")
        break
      }
    }
  }

  # Last resort: scan cache_dir for directories matching index_id
  if (is.null(index_dir)) {
    all_subdirs <- list.dirs(cache_dir, recursive = FALSE, full.names = TRUE)
    for (cand in all_subdirs) {
      if (!grepl(index_id, basename(cand), fixed = TRUE)) next
      bin_files <- list.files(cand, pattern = "\\.bin$", recursive = TRUE)
      if (length(bin_files) > 0) {
        index_dir <- cand
        if (!quiet) message("Found index files at non-standard location: ", cand)
        break
      }
    }
  }

  if (is.null(index_dir)) {
    if (!quiet) {
      message("CRISPRitz index-genome completed but indexed output not detected.")
      message("Expected TST index at: ", tst_dir)
      cache_contents <- list.dirs(cache_dir, recursive = FALSE, full.names = FALSE)
      message("Cache directory contents: ", paste(cache_contents, collapse = ", "))
    }
    warning("Could not locate CRISPRitz indexed genome output.\n",
            "Search will proceed using raw FASTA (slower for bulge mode).")
    return(list(index_dir = genome_dir, indexed = FALSE))
  }

  # ---- 11. Report disk usage ----
  if (!quiet) {
    idx_files_all <- list.files(index_dir, full.names = TRUE, recursive = TRUE)
    idx_size_gb   <- round(sum(file.size(idx_files_all)) / 1e9, 2)
    message("CRISPRitz indexed genome: ", index_dir,
            " (", idx_size_gb, " GB)")
  }

  # ---- 12. Write sentinel ----
  # Store the TST subdirectory path — CRISPRitz parses its basename
  write_unix_lines(c(index_id, index_dir, as.character(Sys.time())), sentinel)

  if (!quiet) message("CRISPRitz index cached successfully.")

  return(list(index_dir = index_dir, indexed = TRUE))
}
