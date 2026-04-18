#' @title CRISPRitz off-target search backend
#'
#' @description Internal functions for off-target searching via the CRISPRitz CLI.
#' Handles input file preparation, CLI invocation (including WSL routing on Windows),
#' and output parsing into a standardised hit data.frame consumed by
#' \code{\link{score_and_aggregate}}.
#'
#' @keywords internal
#' @name ot_crispritz
NULL

#' Locate the CRISPRitz executable
#'
#' Searches for \code{crispritz.py} in the mutateR conda environment, then
#' the system PATH (Unix/macOS), or via WSL (Windows). Used by both
#' \code{\link{run_crispritz_search}} and \code{\link{ensure_crispritz_index}}.
#'
#' @param envname Character. Conda environment name (default \code{"r-mutater"}).
#' @param quiet Logical. Suppress messages (default FALSE).
#'
#' @return A named list:
#'   \describe{
#'     \item{bin}{Character or NULL. Path to the CRISPRitz executable.}
#'     \item{use_wsl}{Logical. Whether execution should be routed through WSL.}
#'     \item{python}{Character or NULL. Path to the conda environment Python
#'           binary (Unix only; NULL when using WSL).}
#'   }
#' @keywords internal
locate_crispritz_bin <- function(envname = "r-mutater", quiet = FALSE) {

  crispritz_bin <- NULL
  use_wsl       <- FALSE
  python_bin    <- NULL

  if (.Platform$OS.type == "unix") {
    # --- Unix/macOS: look in conda env, then PATH ---
    envs <- reticulate::conda_list()
    env_row <- envs[envs$name == envname, ]

    if (nrow(env_row) > 0) {
      env_root  <- dirname(dirname(env_row$python[1]))
      candidate <- file.path(env_root, "bin", "crispritz.py")
      if (file.exists(candidate)) {
        crispritz_bin <- candidate
        python_bin    <- env_row$python[1]
      }
    }

    if (is.null(crispritz_bin)) {
      path_hit <- Sys.which("crispritz.py")
      if (nzchar(path_hit)) crispritz_bin <- path_hit
    }

    if (is.null(python_bin)) python_bin <- Sys.which("python")

  } else {
    # --- Windows: route through WSL ---
    wsl_check <- tryCatch(
      system2("wsl", args = c("bash", "-lc", shQuote("echo ok")),
              stdout = TRUE, stderr = TRUE),
      error = function(e) NULL
    )

    if (!is.null(wsl_check) && any(grepl("ok", wsl_check))) {
      wsl_find <- tryCatch(
        system2("wsl", args = c("bash", "-lc",
                                shQuote("source ~/miniconda3/etc/profile.d/conda.sh && conda activate r-mutater-wsl && which crispritz.py")),
                stdout = TRUE, stderr = TRUE),
        error = function(e) NULL
      )

      if (!is.null(wsl_find) && length(wsl_find) > 0 &&
          nzchar(wsl_find[1]) && !grepl("error|not found", wsl_find[1], ignore.case = TRUE)) {
        crispritz_bin <- trimws(wsl_find[1])
        use_wsl       <- TRUE
      }
    }
  }

  list(bin = crispritz_bin, use_wsl = use_wsl, python = python_bin)
}

#' Copy CRISPRitz index to WSL-native filesystem for faster search
#'
#' WSL accesses Windows NTFS via the 9P protocol, which is 10-100x slower
#' than native ext4 for random-access I/O. CRISPRitz index files (.bin)
#' require intensive random reads. This helper performs a one-time copy
#' of the index directory to the WSL-native filesystem, dramatically
#' accelerating subsequent searches.
#'
#' Results are cached persistently in \code{~/.mutateR_cache/} inside WSL.
#' A sentinel file tracks cache validity. Subsequent calls return the
#' cached native path instantly.
#'
#' Only called on Windows when WSL is in use. No-op on Unix/macOS.
#'
#' @param win_index_dir Character. Windows path to the CRISPRitz TST index
#'   directory (e.g. \code{C:/.../genome_library/NGG_2_...}).
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return Character. WSL-native path to the index directory. On failure,
#'   falls back to the NTFS-mounted WSL path (\code{/mnt/c/...}).
#' @keywords internal
ensure_wsl_native_index <- function(win_index_dir, quiet = FALSE) {

  index_basename <- basename(normalizePath(win_index_dir, winslash = "/",
                                           mustWork = FALSE))

  # ---- 1. Determine WSL home directory ----
  wsl_home <- tryCatch({
    result <- system2("wsl", args = c("bash", "-lc", shQuote("echo $HOME")),
                      stdout = TRUE, stderr = TRUE)
    # Filter to last line starting with / (skip shell startup noise)
    candidates <- grep("^/", trimws(result), value = TRUE)
    if (length(candidates) > 0) candidates[length(candidates)] else NULL
  }, error = function(e) NULL)

  if (is.null(wsl_home) || !nzchar(wsl_home)) {
    if (!quiet) message("Could not determine WSL home directory. Using NTFS path.")
    return(wsl_path(win_index_dir))
  }

  wsl_cache_dir  <- paste0(wsl_home, "/.mutateR_cache")
  wsl_native_dir <- paste0(wsl_cache_dir, "/", index_basename)
  wsl_sentinel   <- paste0(wsl_cache_dir, "/", index_basename, ".native_copy.complete")

  # ---- 2a. Check for existing native copy ----
  check_cmd <- paste0("test -f ", wsl_sentinel, " && echo EXISTS || echo MISSING")
  check_result <- tryCatch(
    system2("wsl", args = c("bash", "-lc", shQuote(check_cmd)),
            stdout = TRUE, stderr = TRUE),
    error = function(e) "MISSING"
  )

  if (any(grepl("EXISTS", check_result))) {
    if (!quiet) message("Using WSL-native index cache: ", wsl_native_dir)
    return(wsl_native_dir)
  }

  # ---- 2b. Reverse migration: old sentinel inside index directory ----
  # Previous versions placed the sentinel inside the index directory.
  # CRISPRitz -index rejects directories containing non-.bin files.
  # Detect the old sentinel, remove it, and write the new external one.
  old_sentinel_path <- paste0(wsl_native_dir, "/.native_copy.complete")
  old_check_cmd <- paste0("test -f ", old_sentinel_path, " && echo OLD_EXISTS || echo OLD_MISSING")
  old_check_result <- tryCatch(
    system2("wsl", args = c("bash", "-lc", shQuote(old_check_cmd)),
            stdout = TRUE, stderr = TRUE),
    error = function(e) "OLD_MISSING"
  )

  if (any(grepl("OLD_EXISTS", old_check_result))) {
    # Verify directory actually contains .bin files before trusting old sentinel
    bin_check_cmd <- paste0("ls ", wsl_native_dir, "/*.bin 2>/dev/null | wc -l")
    bin_count <- tryCatch({
      res <- system2("wsl", args = c("bash", "-lc", shQuote(bin_check_cmd)),
                     stdout = TRUE, stderr = TRUE)
      as.integer(trimws(res[length(res)]))
    }, error = function(e) 0L)

    if (!is.na(bin_count) && bin_count > 0) {
      if (!quiet) message("Migrating WSL-native index sentinel out of index directory ",
                          "(", bin_count, " .bin files verified)...")
      migrate_cmd <- paste0("rm -f ", old_sentinel_path,
                            " && echo COMPLETE > ", wsl_sentinel)
      migrate_result <- tryCatch(
        system2("wsl", args = c("bash", "-lc", shQuote(migrate_cmd)),
                stdout = TRUE, stderr = TRUE),
        error = function(e) { attr(e, "status") <- 1L; e }
      )
      migrate_exit <- attr(migrate_result, "status")
      if (is.null(migrate_exit)) migrate_exit <- 0L

      if (migrate_exit == 0L) {
        if (!quiet) message("Sentinel migrated to: ", wsl_sentinel)
        return(wsl_native_dir)
      } else {
        if (!quiet) message("Sentinel migration failed. Will re-copy index.")
      }
    } else {
      if (!quiet) message("Old sentinel found but no .bin files present. Will re-copy index.")
    }
  }

  # ---- 3. Copy from NTFS to native ext4 ----
  wsl_source <- wsl_path(win_index_dir)

  if (!quiet) {
    # Query source size for progress reporting
    size_cmd <- paste0("du -sh ", wsl_source, " 2>/dev/null | cut -f1")
    size_result <- tryCatch(
      system2("wsl", args = c("bash", "-lc", shQuote(size_cmd)),
              stdout = TRUE, stderr = TRUE),
      error = function(e) NULL
    )
    size_str <- if (!is.null(size_result) && length(size_result) > 0)
      trimws(size_result[1]) else "unknown size"
    message("Copying CRISPRitz index to WSL-native filesystem (", size_str, ")...")
    message("  Source (NTFS): ", wsl_source)
    message("  Target (ext4): ", wsl_native_dir)
    message("  This is a one-time operation (~1-2 min for hg38 index).")
  }

  t_start <- Sys.time()

  copy_cmd <- paste0(
    "mkdir -p ", wsl_cache_dir,
    " && rm -rf ", wsl_native_dir,
    " && cp -r ", wsl_source, " ", wsl_native_dir,
    " && echo COMPLETE > ", wsl_sentinel
  )

  copy_result <- tryCatch(
    system2("wsl", args = c("bash", "-lc", shQuote(copy_cmd)),
            stdout = TRUE, stderr = TRUE),
    error = function(e) {
      attr(e, "status") <- 1L
      e
    }
  )

  exit_code <- attr(copy_result, "status")
  if (is.null(exit_code)) exit_code <- 0L

  t_elapsed <- difftime(Sys.time(), t_start, units = "secs")

  if (exit_code != 0) {
    if (!quiet) {
      message("Failed to copy index to WSL-native filesystem (exit code ", exit_code, ").")
      message("Falling back to NTFS path (significantly slower for CRISPRitz).")
      message("You can manually copy via WSL terminal:")
      message("  mkdir -p ", wsl_cache_dir)
      message("  cp -r ", wsl_source, " ", wsl_native_dir)
    }
    return(wsl_path(win_index_dir))
  }

  if (!quiet) message("Index copied to WSL-native filesystem in ",
                      round(as.numeric(t_elapsed), 0), "s")

  return(wsl_native_dir)
}


#' Reclaim WSL2 page cache between CRISPRitz batches
#'
#' Forces the Linux kernel inside the WSL2 Hyper-V VM to release cached
#' filesystem pages by writing to \code{/proc/sys/vm/drop_caches}. This
#' is necessary because WSL2's VM balloon driver does not proactively
#' release memory back to the host OS, causing the OOM killer to fire
#' prematurely when CRISPRitz burst-loads TST tree files while stale
#' page cache from previous batches occupies the VM's memory ceiling.
#'
#' On native Linux this is a no-op: the kernel reclaims page cache
#' transparently on demand during new allocations, and the OOM killer
#' only fires when genuinely no reclaimable memory remains.
#'
#' The function first attempts \code{drop_caches} without elevated
#' privileges (WSL's default user may be root), then retries with
#' \code{sudo} (many WSL installations have passwordless sudo). Falls
#' back gracefully if neither succeeds.
#'
#' Typically recovers 5--7 GB in the WSL2 VM after several CRISPRitz
#' batches (71 \eqn{\times}{x} ~100 MB \code{.bin} reads per batch).
#' Cost: ~1--2 seconds per invocation.
#'
#' @param quiet Logical. Suppress progress messages (default \code{FALSE}).
#'
#' @return Logical. \code{TRUE} if page cache was successfully dropped;
#'   \code{FALSE} if the operation failed (insufficient permissions or
#'   WSL unavailable). Callers should proceed regardless — failure to
#'   reclaim cache is non-fatal; it simply increases OOM risk for the
#'   next CRISPRitz invocation.
#'
#' @seealso \code{\link{run_crispritz_batched}} (sole caller),
#'   \code{\link{ensure_wsl_native_index}} for related WSL filesystem
#'   helpers.
#'
#' @keywords internal
reclaim_wsl_memory <- function(quiet = FALSE) {

  # ---- 1. Guard: only relevant on Windows/WSL ----
  if (.Platform$OS.type != "windows") return(TRUE)

  # ---- 2. Flush dirty pages and drop page cache via wsl -u root ----
  # wsl -u root runs the command as root inside WSL without requiring
  # sudo configuration. This avoids interactive password prompts and
  # works on any WSL installation regardless of sudoers setup.
  drop_cmd <- "sync && echo 3 > /proc/sys/vm/drop_caches 2>/dev/null"

  result <- tryCatch({
    res <- system2("wsl",
                   args = c("-u", "root", "bash", "-lc", shQuote(drop_cmd)),
                   stdout = TRUE, stderr = TRUE)
    exit <- attr(res, "status")
    if (is.null(exit)) exit <- 0L
    exit == 0L
  }, error = function(e) FALSE)

  if (result) {
    if (!quiet) message("  Reclaimed WSL page cache (drop_caches).")
    return(TRUE)
  }

  # ---- 3. Graceful fallback ----
  if (!quiet) message("  Could not reclaim WSL page cache. ",
                      "Consider adding to .wslconfig:\n",
                      "    [wsl2]\n",
                      "    memory=24GB")
  return(FALSE)
}


#' Write CRISPRitz input files (guides + PAM specification)
#'
#' @param unique_seqs Character vector. Unique protospacer sequences.
#' @param nuclease Character. One of "Cas9", "Cas12a", "enCas12a".
#' @param tmpdir Character. Path to temp directory for file output.
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{guides_file}{Path to the guides text file.}
#'     \item{pam_file}{Path to the PAM specification file.}
#'     \item{n_guides}{Integer. Number of valid guides written.}
#'     \item{guide_length}{Integer. Protospacer length.}
#'     \item{pam_pattern}{Character. Full PAM+guide pattern string.}
#'     \item{clean_seqs}{Character vector. Sanitised protospacers (no PAM) for guide_map key consistency.}
#'     \item{pam_length}{Integer. PAM length (needed by parser to strip PAM from output).}
#'     \item{pam_side}{Character. "3prime" or "5prime".}
#'   }
#' @keywords internal
write_crispritz_inputs <- function(unique_seqs, nuclease, tmpdir) {

  # ---- 1. Validate inputs ----
  if (!is.character(unique_seqs) || length(unique_seqs) == 0) {
    stop("unique_seqs must be a non-empty character vector of protospacer sequences.")
  }

  nuclease <- match.arg(nuclease, c("Cas9", "Cas12a", "enCas12a"))

  if (!dir.exists(tmpdir)) {
    dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- 2. Write guides file ----
  # CRISPRitz expects one protospacer per line, uppercase, no header
  guides_file <- file.path(tmpdir, "guides.txt")

  # Sanitise: uppercase, strip whitespace, validate DNA alphabet
  clean_seqs <- toupper(trimws(unique_seqs))

  invalid <- grep("[^ACGTN]", clean_seqs)
  if (length(invalid) > 0) {
    warning("Removing ", length(invalid),
            " protospacer(s) with non-standard characters: ",
            paste(head(clean_seqs[invalid], 3), collapse = ", "),
            if (length(invalid) > 3) "..." else "")
    clean_seqs <- clean_seqs[-invalid]
  }

  if (length(clean_seqs) == 0) {
    stop("No valid protospacer sequences remaining after filtering.")
  }

  # CRISPRitz requires all guides to be the same length
  seq_lengths <- nchar(clean_seqs)
  if (length(unique(seq_lengths)) > 1) {
    stop("CRISPRitz requires all protospacers to be the same length. ",
         "Found lengths: ", paste(sort(unique(seq_lengths)), collapse = ", "))
  }

  # Pad guides with PAM placeholder to match CRISPRitz's expected input length.
  # CRISPRitz requires guide lines to be the same length as the PAM pattern
  # (spacer positions + PAM positions). Placeholder N's fill the PAM slots.
  pam_placeholder <- switch(
    nuclease,
    "Cas9"     = paste(rep("N", 3), collapse = ""),   # 3' PAM -> append
    "Cas12a"   = paste(rep("N", 4), collapse = ""),   # 5' PAM -> prepend
    "enCas12a" = paste(rep("N", 4), collapse = "")    # 5' PAM -> prepend
  )

  padded_seqs <- switch(
    nuclease,
    "Cas9"     = paste0(clean_seqs, pam_placeholder),
    "Cas12a"   = paste0(pam_placeholder, clean_seqs),
    "enCas12a" = paste0(pam_placeholder, clean_seqs)
  )

  write_unix_lines(padded_seqs, guides_file)

  # ---- 3. Write PAM file ----
  # CRISPRitz PAM format (v2.7.0):
  #   - 3' PAM (Cas9):   sequence is appended after the guide
  #                       e.g. "NGG" with N's in guide positions
  #   - 5' PAM (Cas12a): sequence is prepended before the guide
  #                       e.g. "TTTV" with N's in guide positions
  #
  # File contains a single line with the full PAM+guide-placeholder pattern:
  #   Cas9 (20bp guide, 3' NGG):     NNNNNNNNNNNNNNNNNNNNNGG
  #   Cas12a (23bp guide, 5' TTTV):  TTTvnnnnnnnnnnnnnnnnnnnnnnn
  #
  # Convention: uppercase = PAM positions (exact match required by CRISPRitz)
  #             lowercase n = guide/protospacer positions (variable)

  pam_file <- file.path(tmpdir, "pam.txt")
  guide_len <- seq_lengths[1]
  guide_placeholder <- paste(rep("N", guide_len), collapse = "")

  pam_len <- switch(nuclease,
                    "Cas9"     = 3L,
                    "Cas12a"   = 4L,
                    "enCas12a" = 4L
  )

  pam_string <- switch(nuclease,
                       "Cas9" = {
                         paste0(guide_placeholder, "NGG")
                       },
                       "Cas12a" = {
                         paste0("TTTV", guide_placeholder)
                       },
                       "enCas12a" = {
                         paste0("TTTN", guide_placeholder)
                       }
  )

  # CRISPRitz 2.6.x format: PAM+guide pattern followed by space and PAM length
  pam_line <- paste(pam_string, pam_len)
  write_unix_lines(pam_line, pam_file)

  # ---- 4. Return file paths and metadata ----
  return(list(
    guides_file  = guides_file,
    pam_file     = pam_file,
    n_guides     = length(clean_seqs),
    guide_length = guide_len,
    pam_pattern  = pam_string,
    clean_seqs   = clean_seqs,
    pam_length   = pam_len,
    pam_side     = if (nuclease == "Cas9") "3prime" else "5prime"
  ))
}


#' Execute CRISPRitz genome-wide off-target search
#'
#' @param genome_dir Character. Path to cached FASTA directory.
#' @param index_dir Character or NULL. Path to a pre-indexed CRISPRitz genome
#'   directory, as returned by \code{\link{ensure_crispritz_index}}. When
#'   provided and valid, used in place of \code{genome_dir} for accelerated
#'   bulge-mode search. When NULL (default), the raw FASTA directory is used.
#' @param input_files Named list from \code{write_crispritz_inputs()}.
#' @param max_mismatches Integer.
#' @param max_dna_bulges Integer.
#' @param max_rna_bulges Integer.
#' @param threads Integer.
#' @param envname Character. Conda environment name (default "r-mutater").
#' @param crispritz_info Named list or NULL. Pre-resolved output from
#'   \code{\link{locate_crispritz_bin}}, containing \code{bin}, \code{use_wsl},
#'   and \code{python}. When provided, skips the WSL probe — essential for
#'   batched mode where WSL may degrade after OOM kills. Default NULL
#'   (probe at call time).
#' @param timeout Integer. Maximum seconds for CRISPRitz execution (default 600).
#' @param quiet Logical.
#'
#' @return Character. Path to the CRISPRitz output file.
#' @keywords internal
run_crispritz_search <- function(genome_dir,
                                 input_files,
                                 index_dir = NULL,
                                 max_mismatches,
                                 max_dna_bulges,
                                 max_rna_bulges,
                                 threads = 4L,
                                 envname = "r-mutater",
                                 crispritz_info = NULL,
                                 timeout = 600L,
                                 quiet = FALSE) {

  # ---- 1. Locate CRISPRitz executable ----
  # Accept pre-resolved CRISPRitz info to avoid re-probing WSL per batch.
  # WSL can enter a degraded state after OOM kills; re-probing then fails
  # with "CRISPRitz not found" even though the installation is intact.
  cz <- if (!is.null(crispritz_info)) crispritz_info else
    locate_crispritz_bin(envname = envname, quiet = quiet)

  if (is.null(cz$bin)) {
    if (.Platform$OS.type == "unix") {
      stop("CRISPRitz executable not found.\n",
           "Looked in conda env '", envname, "' and system PATH.\n",
           "Please run install_mutater_env(fresh = TRUE) to install CRISPRitz.")
    } else {
      stop("CRISPRitz not found via WSL.\n",
           "On Windows, off-target scoring requires CRISPRitz installed in WSL.\n",
           "Setup instructions:\n",
           "  1. Open a WSL terminal\n",
           "  2. conda create -n r-mutater-wsl python=3.9\n",
           "  3. conda install -n r-mutater-wsl -c bioconda -c conda-forge crispritz\n",
           "Then re-run the mutateR pipeline.")
    }
  }

  crispritz_bin <- cz$bin
  use_wsl       <- cz$use_wsl

  if (!quiet) message("Using CRISPRitz: ", crispritz_bin,
                      if (use_wsl) " (via WSL)" else "")

  # ---- 2a. Determine search mode ----
  use_bulges <- (max_dna_bulges > 0 || max_rna_bulges > 0)

  # ---- 2b. Resolve effective genome directory ----
  # Use pre-indexed directory if available; fall back to raw FASTA.
  # CRISPRitz search -index expects the TST index directory directly —
  # it parses the directory basename to extract PAM and max bulge count
  # (convention: {PAM}_{maxBulges}_{genomeName}).
  effective_genome_dir <- genome_dir
  use_index <- FALSE

  if (!is.null(index_dir) && dir.exists(index_dir)) {
    effective_genome_dir <- index_dir
    use_index <- TRUE
    if (!quiet) message("Using pre-indexed genome: ", index_dir)
  }

  # ---- 3. Construct output prefix ----
  output_prefix <- file.path(dirname(input_files$guides_file), "crispritz_results")

  # ---- 4. Build command arguments ----
  # CRISPRitz 2.6.x CLI:
  #   crispritz.py search <genome_dir> <pam_file> <guides_file> <output>
  #     -mm <N> [-bDNA <N> -bRNA <N>] [-th <N>] -r

  args <- c(
    "search",
    effective_genome_dir,
    input_files$pam_file,
    input_files$guides_file,
    output_prefix,
    "-mm", as.character(max_mismatches)
  )

  if (use_bulges) {
    args <- c(args,
              "-bDNA", as.character(max_dna_bulges),
              "-bRNA", as.character(max_rna_bulges))
  }

  if (use_index) args <- c(args, "-index")
  args <- c(args, "-th", as.character(threads), "-r")

  if (!quiet) {
    message("CRISPRitz mode: ",
            if (use_index) "indexed search (pre-built index)" else "de novo search (no pre-built index)")
    if (!use_wsl) {
      message("CRISPRitz command: crispritz.py ", paste(args, collapse = " "))
    }
    message("Searching ", input_files$n_guides,
            " unique guides (mm=", max_mismatches,
            if (use_bulges) paste0(", bDNA=", max_dna_bulges,
                                   ", bRNA=", max_rna_bulges) else "",
            ")...")
  }

  # ---- 5. Execute ----
  stdout_log <- file.path(dirname(output_prefix), "crispritz_stdout.log")
  stderr_log <- file.path(dirname(output_prefix), "crispritz_stderr.log")

  if (use_wsl) {
    if (use_index) {
      wsl_genome_dir <- ensure_wsl_native_index(effective_genome_dir, quiet = quiet)
    } else {
      wsl_genome_dir <- wsl_path(effective_genome_dir)
    }
    wsl_pam_file   <- wsl_path(input_files$pam_file)
    wsl_guides     <- wsl_path(input_files$guides_file)
    wsl_output     <- wsl_path(output_prefix)

    wsl_args <- c(
      "search", wsl_genome_dir, wsl_pam_file, wsl_guides, wsl_output,
      "-mm", as.character(max_mismatches)
    )
    if (use_bulges) {
      wsl_args <- c(wsl_args,
                    "-bDNA", as.character(max_dna_bulges),
                    "-bRNA", as.character(max_rna_bulges))
    }

    if (use_index) wsl_args <- c(wsl_args, "-index")
    wsl_args <- c(wsl_args, "-th", as.character(threads), "-r")

    wsl_stdout <- wsl_path(stdout_log)
    wsl_stderr <- wsl_path(stderr_log)

    wsl_cmd_str <- paste(
      "source ~/miniconda3/etc/profile.d/conda.sh",
      "&&", "conda activate r-mutater-wsl",
      "&&", "timeout", as.character(timeout),
      "crispritz.py", paste(wsl_args, collapse = " "),
      ">", wsl_stdout,
      "2>", wsl_stderr
    )

    if (!quiet) message("CRISPRitz command (WSL): crispritz.py ",
                        paste(wsl_args, collapse = " "))

    exit_code <- system2(
      "wsl",
      args = c("bash", "-lc", shQuote(wsl_cmd_str)),
      stdout = TRUE,
      stderr = TRUE
    )

    exit_code <- attr(exit_code, "status")
    if (is.null(exit_code)) exit_code <- 0L

  } else {
    # Unix/macOS: direct execution
    if (grepl("\\.py$", crispritz_bin)) {
      exec <- cz$python
      full_args <- c(crispritz_bin, args)
    } else {
      exec <- crispritz_bin
      full_args <- args
    }

    exit_code <- system2(
      command = exec,
      args = full_args,
      stdout = stdout_log,
      stderr = stderr_log
    )
  }

  # ---- 5b. Diagnostic logging (non-quiet mode) ----
  if (!quiet) {
    stderr_content <- if (file.exists(stderr_log)) {
      tryCatch(readLines(stderr_log, warn = FALSE, n = 50),
               error = function(e) "(could not read stderr log)")
    } else "(stderr log not created)"

    stdout_content <- if (file.exists(stdout_log)) {
      tryCatch(readLines(stdout_log, warn = FALSE, n = 50),
               error = function(e) "(could not read stdout log)")
    } else "(stdout log not created)"

    if (length(stdout_content) > 0 && any(nzchar(stdout_content))) {
      message("CRISPRitz stdout:\n  ", paste(stdout_content, collapse = "\n  "))
    }
    if (length(stderr_content) > 0 && any(nzchar(stderr_content))) {
      message("CRISPRitz stderr:\n  ", paste(stderr_content, collapse = "\n  "))
    }

    out_dir_contents <- list.files(dirname(output_prefix), full.names = FALSE)
    message("Output directory contents: ", paste(out_dir_contents, collapse = ", "))
  }

  # ---- 5c. Detect CRISPRitz errors reported via stdout (exit code may still be 0) ----
  timeout_fallback <- FALSE

  stdout_lines <- if (file.exists(stdout_log)) {
    tryCatch(readLines(stdout_log, warn = FALSE),
             error = function(e) character(0))
  } else character(0)

  if (any(grepl("^ERROR", stdout_lines, ignore.case = TRUE))) {
    err_detail <- paste(grep("ERROR", stdout_lines, ignore.case = TRUE, value = TRUE),
                        collapse = "\n")
    stop("CRISPRitz reported an error:\n", err_detail)
  }

  # ---- 5d. Handle timeout with partial-results fallback ----
  if (exit_code == 124) {
    # CRISPRitz may have completed the mismatch pass (and partial bulge passes)
    # before timing out during the merge step. Check for usable partial output.
    partial_candidates <- list.files(
      dirname(output_prefix),
      pattern = paste0("^", basename(output_prefix),
                       ".*\\.(targets|bestMerge|altMerge)\\.txt$"),
      full.names = TRUE
    )
    partial_with_content <- partial_candidates[
      vapply(partial_candidates, function(f) file.size(f) > 0, logical(1))
    ]

    if (length(partial_with_content) > 0) {
      n_lines <- length(readLines(partial_with_content[1], warn = FALSE))
      warning("CRISPRitz search timed out after ", timeout, " seconds, ",
              "but partial results are available (", n_lines, " hits).\n",
              "Using partial output. Note:\n",
              "  - Results may contain duplicate hits for the same locus ",
              "(conservative: overestimates off-target burden).\n",
              "  - Some guides/chromosomes may be incompletely searched.\n",
              "Consider increasing ot_timeout or using ot_backend='hybrid'.")
      timeout_fallback <- TRUE
    } else {
      stop("CRISPRitz search timed out after ", timeout, " seconds ",
           "with no usable output.\n",
           "Search mode: ",
           if (use_index) "indexed" else "de novo (no pre-built index)", "\n",
           if (!use_index) {
             paste0("A pre-built index may dramatically reduce search time.\n",
                    "Run ensure_crispritz_index() separately, ",
                    "or use ot_backend='hybrid'.\n")
           } else "",
           "Consider reducing max_mismatches or the number of input gRNAs.")
    }
  }

  # ---- 6. Check execution status ----
  if (exit_code != 0 && !timeout_fallback) {
    err_msg <- tryCatch(
      paste(readLines(stderr_log, warn = FALSE), collapse = "\n"),
      error = function(e) "(could not read stderr log)"
    )

    stop("CRISPRitz search failed (exit code ", exit_code, ").\n",
         "Search mode: ", if (use_index) "indexed" else "de novo", "\n",
         "Genome dir: ", effective_genome_dir, "\n",
         "Command: crispritz.py ", paste(args, collapse = " "), "\n",
         "Stderr:\n", err_msg, "\n",
         "If CRISPRitz was found on system PATH rather than in the conda env,\n",
         "ensure all CRISPRitz dependencies are available.\n",
         "Recommended: install_mutater_env(fresh = TRUE)")
  }

  if (!quiet) {
    if (timeout_fallback) {
      message("CRISPRitz search: using partial results after timeout.")
    } else {
      message("CRISPRitz search completed successfully.")
    }
  }

  # ---- 7. Locate output file(s) ----
  # CRISPRitz writes results with predictable naming:
  #   {output_prefix}.targets.txt         (mismatch-only mode)
  #   {output_prefix}.bestMerge.txt       (bulge mode — merged best hits)
  #   {output_prefix}.altMerge.txt        (bulge mode — alternative hits)

  if (use_bulges) {
    result_file <- paste0(output_prefix, ".bestMerge.txt")
  } else {
    result_file <- paste0(output_prefix, ".targets.txt")
  }

  # Fallback: search for any result file matching the prefix
  if (!file.exists(result_file)) {
    candidates <- list.files(
      dirname(output_prefix),
      pattern = paste0("^", basename(output_prefix), ".*\\.(targets|bestMerge|altMerge)\\.txt$"),
      full.names = TRUE
    )

    if (length(candidates) == 0) {
      stop("CRISPRitz completed but no output file found.\n",
           "Expected: ", result_file, "\n",
           "Directory contents: ",
           paste(list.files(dirname(output_prefix)), collapse = ", "))
    }

    priority <- c("bestMerge", "targets", "altMerge")
    for (p in priority) {
      matched_files <- grep(p, candidates, value = TRUE)
      if (length(matched_files) > 0) {
        result_file <- matched_files[1]
        break
      }
    }

    if (!quiet) message("Using output file: ", basename(result_file))
  }

  # ---- 8a. Quick sanity check on output ----
  if (file.size(result_file) == 0) {
    warning("CRISPRitz output file is empty — no off-target hits found. ",
            "All gRNAs will receive maximum specificity scores.")
  }

  if (timeout_fallback) {
    attr(result_file, "partial") <- TRUE
    attr(result_file, "timeout")  <- timeout
  }

  # ---- 8b. Output size guard (prevents multi-GB dedup) ----
  max_output_bytes <- getOption("mutateR.max_crispritz_output_bytes",
                                default = 200 * 1024 * 1024)  # 200 MB default
  raw_size <- file.size(result_file)

  if (!is.na(raw_size) && raw_size > max_output_bytes) {
    stop("CRISPRitz output exceeds size limit (",
         round(raw_size / 1024 / 1024), " MB > ",
         round(max_output_bytes / 1024 / 1024), " MB). ",
         "Guide likely targets repetitive sequence.")
  }

  # ---- 8c. Disk-based deduplication for indexed-mode output ----
  # CRISPRitz -index mode writes all hit representations (mismatch, DNA bulge,
  # RNA bulge) for the same locus without merging. This creates
  # redundancy that inflates R's peak memory during parsing. Dedup on disk
  # via awk before R reads the file, reducing peak parse memory.
  #
  # Non-indexed mode (.bestMerge.txt or brute-force .targets.txt) is already
  # deduplicated by CRISPRitz's internal merge step — skip.
  # R-side dedup in parse_crispritz_output() step 8b is retained as safety net.
  if (use_index) {
    result_file <- dedup_crispritz_on_disk(
      result_file = result_file,
      use_wsl     = use_wsl,
      quiet       = quiet
    )
  }

  return(result_file)
}


#' Generate awk deduplication script for CRISPRitz indexed output
#'
#' Returns the awk dedup script as a character vector (one element per line).
#' The script deduplicates CRISPRitz \code{.targets.txt} output from
#' \code{-index} mode, where the same (guide, chr, pos, strand) locus
#' appears multiple times with different alignment representations
#' (mismatch, DNA bulge, RNA bulge).
#'
#' For each unique locus, keeps the representation with the lowest total
#' edits, then lowest bulge size, then prefers mismatch-only (type "X").
#' This mirrors the R-side dedup logic in \code{parse_crispritz_output()}
#' but operates on disk before R reads the file.
#'
#' The script is generated at runtime rather than shipped as an
#' \code{inst/} file because R's package build toolchain and custom I/O
#' helpers (e.g. \code{write_unix_lines()}) corrupt awk's
#' \code{$}-column references. Standard R string literals handle
#' \code{$} without issue; \code{writeLines()} on a binary connection
#' writes them faithfully.
#'
#' @return Character vector. One element per line of the awk script.
#'
#' @details
#' CRISPRitz \code{.targets.txt} columns (tab-delimited, no header):
#' \enumerate{
#'   \item bulge_type ("X", "DNA", or "RNA")
#'   \item crRNA (may contain gap characters "-")
#'   \item DNA (off-target genomic sequence, may contain gaps)
#'   \item chromosome
#'   \item position (0-based)
#'   \item cluster_position
#'   \item direction ("+" or "-")
#'   \item mismatches
#'   \item bulge_size
#'   \item total (mismatches + bulge_size)
#'   \item PAM_gen
#'   \item guide_with_PAM
#' }
#'
#' Dedup key: gap-stripped crRNA (\code{\$2}), chromosome (\code{\$4}),
#' position (\code{\$5}), strand (\code{\$7}). Joined via awk's
#' \code{SUBSEP} (non-printable, no collision risk).
#'
#' Preference (lower wins): total_edits (\code{\$10}), then bulge_size
#' (\code{\$9}), then bulge_type rank (X=0, DNA=1, RNA=2, other=3).
#'
#' Stats line written to stderr: \code{<input_rows> <output_rows>}
#' (space-separated integers), parsed by \code{dedup_crispritz_on_disk()}.
#'
#' @keywords internal
generate_dedup_awk <- function() {

  # ---- Placeholder strategy ----
  # Awk field references (\$1, \$2, etc.) are valid in R string literals,

  # but intermediaries (e.g. clipboard, write_unix_lines, R CMD INSTALL)
  # routinely corrupt $ → \$ during transmission. \$ is not a valid R
  # escape sequence, causing parse errors at source time.
  #
  # Solution: use @ as a placeholder for $ throughout. @ has no meaning
  # in awk syntax. A single gsub("@", "$", ..., fixed = TRUE) at the
  # end produces the correct awk code. This is immune to any
  # intermediary that escapes $.

  lines <- c(
    "BEGIN {",
    "    FS = \"\\t\"",
    "    OFS = \"\\t\"",
    "}",

    "NF < 10 { print; next }",

    "{",
    "    guide = @2",
    "    gsub(\"-\", \"\", guide)",
    "    key = guide SUBSEP @4 SUBSEP @5 SUBSEP @7",
    "",
    "    te = @10 + 0",
    "    bs = @9 + 0",
    "    if (@1 == \"X\") bt = 0",
    "    else if (@1 == \"DNA\") bt = 1",
    "    else if (@1 == \"RNA\") bt = 2",
    "    else bt = 3",
    "",
    "    if (!(key in best_te)) {",
    "        best_te[key] = te",
    "        best_bs[key] = bs",
    "        best_bt[key] = bt",
    "        best_line[key] = @0",
    "    } else {",
    "        if (te < best_te[key] ||",
    "            (te == best_te[key] && bs < best_bs[key]) ||",
    "            (te == best_te[key] && bs == best_bs[key] && bt < best_bt[key])) {",
    "            best_te[key] = te",
    "            best_bs[key] = bs",
    "            best_bt[key] = bt",
    "            best_line[key] = @0",
    "        }",
    "    }",
    "}",

    "END {",
    "    n = 0",
    "    for (key in best_line) {",
    "        print best_line[key]",
    "        n++",
    "    }",
    "    printf \"%d %d\\n\", NR, n > \"/dev/stderr\"",
    "}"
  )

  gsub("@", "$", lines, fixed = TRUE)
}

#' Write an awk script to disk with verified integrity
#'
#' Writes a character vector of awk code to a file with guaranteed Unix
#' line endings (LF) and no character transformation. After writing,
#' reads back the file and verifies that awk field references (\code{$})
#' survived intact — a defensive check against the intermediary corruption
#' that has plagued every prior delivery mechanism.
#'
#' Uses a binary connection to prevent Windows \code{\\r\\n} translation.
#' Does not use \code{write_unix_lines()} (which escapes \code{$}).
#'
#' @param lines Character vector. One element per line of awk code,
#'   as returned by \code{\link{generate_dedup_awk}}.
#' @param path Character. Output file path.
#'
#' @return Logical. \code{TRUE} if the file was written and verified
#'   successfully; \code{FALSE} on any failure (write error, missing
#'   \code{$} references, etc.). On \code{FALSE}, the caller should
#'   fall back to R-side deduplication.
#'
#' @keywords internal
write_awk_script <- function(lines, path) {

  # ---- 1. Write via binary connection (LF line endings, no transformation) ----
  ok <- tryCatch({
    con <- file(path, open = "wb")
    writeLines(lines, con, sep = "\n")
    close(con)
    TRUE
  }, error = function(e) {
    tryCatch(close(con), error = function(e2) NULL)
    warning("Failed to write awk script to ", path, ": ", e$message)
    FALSE
  })

  if (!ok) return(FALSE)

  # ---- 2. Verify file exists and is non-empty ----
  if (!file.exists(path) || file.size(path) == 0) {
    warning("Awk script file is missing or empty after write: ", path)
    return(FALSE)
  }

  # ---- 3. Read back and verify $ references survived ----
  readback <- tryCatch(
    readLines(path, n = 20, warn = FALSE),
    error = function(e) NULL
  )

  if (is.null(readback) || length(readback) == 0) {
    warning("Could not read back awk script for verification: ", path)
    return(FALSE)
  }

  has_dollar <- any(grepl("$", readback, fixed = TRUE))
  if (!has_dollar) {
    warning("Awk script verification failed: no $ field references found. ",
            "File may be corrupted. Contents:\n",
            paste(head(readback, 10), collapse = "\n"))
    return(FALSE)
  }

  # ---- 4. Spot-check: a known line should match exactly ----
  has_clean_ref <- any(grepl("guide = $", readback, fixed = TRUE))
  if (!has_clean_ref) {
    warning("Awk script verification failed: expected 'guide = $2' but ",
            "pattern not found. $ may have been escaped. First 10 lines:\n",
            paste(head(readback, 10), collapse = "\n"))
    return(FALSE)
  }

  return(TRUE)
}


#' Deduplicate CRISPRitz indexed output on disk before R parsing
#'
#' CRISPRitz \code{-index} mode writes all hit representations (mismatch, DNA
#' bulge, RNA bulge) for the same genomic locus to \code{.targets.txt} without
#' the merge step that normally produces \code{.bestMerge.txt}. This creates
#' 70--82\% redundancy: the same (guide, chr, pos, strand) tuple appears
#' multiple times with different alignment types.
#'
#' This function runs a single-pass \code{awk} script on the raw output file
#' to deduplicate hits on disk \emph{before} R reads them, reducing peak parse
#' memory.
#'
#' Only called for indexed-mode \code{.targets.txt} output. Non-indexed
#' \code{.targets.txt} and \code{.bestMerge.txt} files are already
#' deduplicated by CRISPRitz's internal merge step.
#'
#' @param result_file Character. Path to CRISPRitz \code{.targets.txt} output.
#' @param use_wsl Logical. Whether \code{awk} should be invoked via WSL
#'   (Windows). On Unix/macOS, \code{awk} is called directly.
#' @param quiet Logical. Suppress progress messages (default FALSE).
#'
#' @return Character. Path to the deduplicated file (sibling of the original
#'   with \code{.dedup} suffix). On failure, returns the original path
#'   unchanged with a warning. Attributes (\code{partial}, \code{timeout})
#'   from the original path are transferred to the deduplicated path.
#'
#' @details The \code{awk} script is written to a temporary file in the same
#'   directory as \code{result_file} to avoid shell quoting issues. It uses
#'   \code{SUBSEP} (a non-printable character) as the composite key separator
#'   and outputs row counts to \code{/dev/stderr} for diagnostic reporting.
#'
#'   The R-side dedup in \code{parse_crispritz_output()} step 8b is retained
#'   as a safety net: if disk dedup fails silently or partially, R catches
#'   any remaining duplicates. When disk dedup succeeds, the R-side pass
#'   finds zero duplicates and exits in O(n).
#'
#' @keywords internal
dedup_crispritz_on_disk <- function(result_file, use_wsl = FALSE, quiet = FALSE) {

  # ---- 1. Guard: skip if file is empty or missing ----
  if (!file.exists(result_file) || file.size(result_file) == 0) {
    return(result_file)
  }

  # ---- 2. Generate awk dedup script at runtime ----
  # Generated in-memory rather than shipped as inst/awk/dedup_indexed.awk
  # because R's package build toolchain and write_unix_lines() corrupt
  # awk's $-column references. generate_dedup_awk() uses a placeholder
  # substitution strategy that is immune to intermediary escaping.
  awk_lines <- generate_dedup_awk()

  # ---- 3. Write awk script to working directory ----
  script_file <- file.path(dirname(result_file), "dedup_indexed.awk")
  script_ok <- write_awk_script(awk_lines, script_file)
  on.exit(unlink(script_file), add = TRUE)

  if (!script_ok) {
    if (!quiet) message("Could not write verified awk dedup script. ",
                        "Falling back to R-side dedup.")
    return(result_file)
  }

  # ---- 4. Set up output path ----
  dedup_file <- paste0(result_file, ".dedup")

  # ---- 5. Execute ----
  # Two-file approach: the awk logic lives in dedup_indexed.awk (from inst/),
  # and a tiny shell wrapper script wires up the I/O paths. Both are written
  # to disk to completely bypass multi-layer quoting (system2 → Windows →
  # WSL → bash -c) that corrupts redirect operators on Windows.
  t_start     <- Sys.time()
  exec_result <- NULL

  if (use_wsl) {
    wsl_script <- wsl_path(script_file)
    wsl_input  <- wsl_path(result_file)
    wsl_output <- wsl_path(dedup_file)

    # Stats file captures awk's stderr (row counts for diagnostics).
    # Placed alongside the result file for automatic tmpdir cleanup.
    stats_file <- file.path(dirname(result_file), "dedup_stats.txt")

    wrapper_file <- file.path(dirname(result_file), "dedup_wrapper.sh")
    write_unix_lines(c(
      "#!/bin/bash",
      paste("awk -f", wsl_script, wsl_input,
            ">", wsl_output,
            "2>", wsl_path(stats_file))
    ), wrapper_file)
    on.exit(unlink(wrapper_file), add = TRUE)
    on.exit(unlink(stats_file), add = TRUE)

    wsl_wrapper <- wsl_path(wrapper_file)

    exec_result <- tryCatch(
      system2("wsl", args = c("bash", wsl_wrapper),
              stdout = TRUE, stderr = TRUE),
      error = function(e) { attr(e, "status") <- 1L; e }
    )
    exit_code <- attr(exec_result, "status")
    if (is.null(exit_code)) exit_code <- 0L

  } else {
    # Direct Unix/macOS execution — no quoting issues
    stats_file <- file.path(dirname(result_file), "dedup_stats.txt")
    on.exit(unlink(stats_file), add = TRUE)

    exit_code <- system2("awk",
                         args = c("-f", script_file, result_file),
                         stdout = dedup_file,
                         stderr = stats_file)
  }

  t_elapsed <- round(difftime(Sys.time(), t_start, units = "secs"), 1)

  # ---- 6. Validate output ----
  if (exit_code != 0 || !file.exists(dedup_file) || file.size(dedup_file) == 0) {
    if (!quiet) {
      message("Disk-based dedup failed (exit code ", exit_code,
              "). Falling back to R-side dedup.")
      # Diagnostic: print captured system2 output (bash/WSL errors)
      if (is.character(exec_result) && length(exec_result) > 0 &&
          any(nzchar(exec_result))) {
        message("  WSL output: ",
                paste(utils::head(exec_result, 10), collapse = "\n  "))
      }
      # Diagnostic: print stats/stderr file (awk errors)
      if (file.exists(stats_file) && file.size(stats_file) > 0) {
        stats_content <- tryCatch(
          readLines(stats_file, n = 5, warn = FALSE),
          error = function(e) NULL
        )
        if (!is.null(stats_content) && length(stats_content) > 0) {
          message("  awk stderr: ",
                  paste(stats_content, collapse = "\n  "))
        }
      }
      # Diagnostic: print the wrapper script for manual reproduction
      if (use_wsl && file.exists(wrapper_file)) {
        wrapper_content <- tryCatch(readLines(wrapper_file, warn = FALSE),
                                    error = function(e) NULL)
        if (!is.null(wrapper_content)) {
          message("  Wrapper script contents:\n    ",
                  paste(wrapper_content, collapse = "\n    "))
          message("  To reproduce manually: wsl bash ", wsl_wrapper)
        }
      }
    }
    unlink(dedup_file)
    return(result_file)
  }

  # ---- 7. Sanity check: dedup output must not be larger than input ----
  size_before <- file.size(result_file)
  size_after  <- file.size(dedup_file)

  if (size_after > size_before) {
    if (!quiet) message("WARNING: Dedup output larger than input. ",
                        "Using original file (R-side dedup will handle).")
    unlink(dedup_file)
    return(result_file)
  }

  # ---- 8. Parse stats and report ----
  n_before <- NA_integer_
  n_after  <- NA_integer_

  if (file.exists(stats_file) && file.size(stats_file) > 0) {
    stats_line <- tryCatch(
      readLines(stats_file, n = 1L, warn = FALSE),
      error = function(e) NULL
    )
    if (!is.null(stats_line) && nzchar(stats_line)) {
      parts <- suppressWarnings(
        as.integer(strsplit(trimws(stats_line), "\\s+")[[1]])
      )
      if (length(parts) >= 2 && !any(is.na(parts))) {
        n_before <- parts[1]
        n_after  <- parts[2]
      }
    }
  }

  if (!quiet) {
    if (!is.na(n_before) && !is.na(n_after)) {
      n_removed     <- n_before - n_after
      reduction_pct <- round(100 * n_removed / n_before, 1)
      message(sprintf(
        "Disk-based dedup: %s -> %s rows (%s removed, %s%%) in %ss",
        format(n_before, big.mark = ","),
        format(n_after, big.mark = ","),
        format(n_removed, big.mark = ","),
        reduction_pct,
        t_elapsed
      ))
    } else {
      message("Disk-based dedup completed in ", t_elapsed, "s ",
              "(row counts unavailable; file size: ",
              round(size_before / 1024, 1), " KB -> ",
              round(size_after / 1024, 1), " KB)")
    }
  }

  # ---- 9. Transfer attributes from original path ----
  # run_crispritz_search() may have set 'partial' and 'timeout' attributes
  # on the original result_file path for timeout-fallback handling.
  # These must be preserved on the deduplicated path so downstream
  # consumers (parse_crispritz_output, score_and_aggregate) see them.
  orig_partial <- attr(result_file, "partial")
  orig_timeout <- attr(result_file, "timeout")
  if (!is.null(orig_partial)) attr(dedup_file, "partial") <- orig_partial
  if (!is.null(orig_timeout)) attr(dedup_file, "timeout") <- orig_timeout

  # ---- 10. Clean up original (large) file to free disk space ----
  unlink(result_file)

  return(dedup_file)
}


#' Parse CRISPRitz output into a structured data frame
#'
#' Reads CRISPRitz tab-delimited output and maps hits back to input gRNAs
#' via the guide_map index. PAM sequences are stripped from both guide and
#' off-target sequences so that downstream CFD scoring operates on
#' protospacer-length strings.
#'
#' @param output_path Character. Path to CRISPRitz results file.
#' @param guide_map Named list. Maps unique protospacer sequence (names) to
#'        vectors of original GRanges indices (values).
#' @param pam_length Integer. Length of the PAM sequence (default 3).
#' @param pam_side Character. "3prime" or "5prime" (default "3prime").
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{guide_seq}{Character. The protospacer sequence that was searched.}
#'     \item{grna_indices}{List-column. Vector of original GRanges indices sharing this protospacer.}
#'     \item{chr}{Character. Chromosome of the off-target hit.}
#'     \item{pos}{Integer. Genomic position (1-based).}
#'     \item{strand}{Character. "+" or "-".}
#'     \item{offtarget_seq}{Character. Genomic off-target sequence (PAM-stripped).}
#'     \item{pam_gen}{Character. PAM sequence found in genome.}
#'     \item{n_mismatches}{Integer. Number of mismatches.}
#'     \item{bulge_type}{Character. "X" (mismatch only), "DNA", or "RNA".}
#'     \item{bulge_size}{Integer. Size of the bulge (0 for mismatch-only).}
#'   }
#'   Returns a zero-row data frame with correct column types if no hits are present.
#' @keywords internal
parse_crispritz_output <- function(output_path, guide_map, pam_length = 3L, pam_side = "3prime") {

  # ---- 1. Define canonical output schema ----
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

  # ---- 2. Handle empty output ----
  if (!file.exists(output_path) || file.size(output_path) == 0) {
    return(empty_result)
  }

  # ---- 3. Read CRISPRitz output ----
  # CRISPRitz 2.7.0 tab-delimited output columns (bestMerge / targets):
  #   1: bulge_type       — "X" (mismatch), "DNA", or "RNA"
  #   2: crRNA            — the guide/spacer sequence (with gaps if bulge)
  #   3: DNA              — the off-target genomic sequence (with gaps if bulge)
  #   4: chromosome       — chromosome name
  #   5: position         — 0-based genomic start position
  #   6: cluster_position — cluster position (used internally by CRISPRitz)
  #   7: direction        — "+" or "-"
  #   8: mismatches       — number of mismatches
  #   9: bulge_size       — size of the bulge (0 for mismatch-only)
  #  10: total            — total edit distance (mismatches + bulge_size)
  #  11: PAM_gen          — PAM sequence found in genome
  #  12: guide_with_PAM   — full guide+PAM or PAM+guide string

  raw <- tryCatch(
    utils::read.delim(
      output_path,
      header = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "#",
      quote = ""
    ),
    error = function(e) {
      warning("Failed to read CRISPRitz output: ", e$message)
      return(NULL)
    }
  )

  if (is.null(raw) || nrow(raw) == 0) {
    return(empty_result)
  }

  # ---- 4. Assign column names based on detected width ----
  expected_cols <- c("bulge_type", "crRNA", "DNA", "chromosome", "position",
                     "cluster_position", "direction", "mismatches",
                     "bulge_size", "total", "PAM_gen", "guide_with_PAM")

  if (ncol(raw) >= length(expected_cols)) {
    names(raw)[seq_along(expected_cols)] <- expected_cols
  } else if (ncol(raw) >= 10) {
    names(raw)[1:10] <- expected_cols[1:10]
  } else {
    warning("CRISPRitz output has unexpected format (", ncol(raw),
            " columns). Expected >= 10.")
    return(empty_result)
  }

  # ---- 5. Clean and type-cast ----
  hits <- data.frame(
    bulge_type    = as.character(raw$bulge_type),
    crRNA_raw     = as.character(raw$crRNA),
    offtarget_seq = as.character(raw$DNA),
    chr           = as.character(raw$chromosome),
    pos           = as.integer(raw$position),
    strand        = as.character(raw$direction),
    n_mismatches  = as.integer(raw$mismatches),
    bulge_size    = as.integer(raw$bulge_size),
    total_edits   = if ("total" %in% names(raw)) as.integer(raw$total) else NA_integer_,
    pam_gen       = if ("PAM_gen" %in% names(raw)) toupper(as.character(raw$PAM_gen)) else NA_character_,
    stringsAsFactors = FALSE
  )

  # Free the raw read.delim() result
  rm(raw); gc(verbose = FALSE)

  # Sanitise strand values — CRISPRitz may output padded or non-standard values
  hits$strand <- trimws(hits$strand)
  hits$strand[!hits$strand %in% c("+", "-")] <- "*"

  # CRISPRitz positions are 0-based; convert to 1-based for R/Bioconductor
  hits$pos <- hits$pos + 1L

  # ---- 6. Extract clean guide sequence for mapping ----
  # Strip gaps, then remove PAM placeholder to recover the original protospacer
  hits$guide_seq_full <- toupper(gsub("-", "", hits$crRNA_raw))
  hits$guide_seq <- if (pam_side == "3prime") {
    substr(hits$guide_seq_full, 1, nchar(hits$guide_seq_full) - pam_length)
  } else {
    substr(hits$guide_seq_full, pam_length + 1, nchar(hits$guide_seq_full))
  }

  # Also strip gaps from off-target sequence for downstream CFD scoring
  hits$offtarget_seq_clean <- toupper(gsub("-", "", hits$offtarget_seq))

  # Then strip PAM from off-target sequence to match protospacer length
  hits$offtarget_seq_clean <- if (pam_side == "3prime") {
    substr(hits$offtarget_seq_clean, 1, nchar(hits$offtarget_seq_clean) - pam_length)
  } else {
    substr(hits$offtarget_seq_clean, pam_length + 1, nchar(hits$offtarget_seq_clean))
  }

  # ---- 7. Map guide sequences to GRanges indices ----
  hits$grna_indices <- lapply(hits$guide_seq, function(seq) {
    idx <- guide_map[[seq]]
    if (is.null(idx)) integer(0) else idx
  })

  unmapped <- vapply(hits$grna_indices, function(x) length(x) == 0, logical(1))
  if (any(unmapped)) {
    n_unmapped <- sum(unmapped)
    warning(n_unmapped, " off-target hit(s) could not be mapped back to input gRNAs. ",
            "These will be excluded from scoring.")
    hits <- hits[!unmapped, , drop = FALSE]
  }

  if (nrow(hits) == 0) {
    return(empty_result)
  }

  # ---- 8a. Standardise bulge_type labels ----
  # CRISPRitz can report combined bulge types for hits with both DNA and RNA
  # bulges, e.g. "RNA,DNA" or "DNA,RNA". We normalise these to "DNA+RNA"
  # (a combined type treated as a bulge hit for scoring purposes).
  hits$bulge_type <- toupper(trimws(hits$bulge_type))

  # Normalise combined types: "RNA,DNA" / "DNA,RNA" -> "DNA+RNA"
  is_combined <- grepl(",", hits$bulge_type, fixed = TRUE)
  if (any(is_combined)) {
    hits$bulge_type[is_combined] <- "DNA+RNA"
  }

  valid_types <- c("X", "DNA", "RNA", "DNA+RNA")
  invalid_type <- !hits$bulge_type %in% valid_types
  if (any(invalid_type)) {
    warning("Unexpected bulge_type values: ",
            paste(unique(hits$bulge_type[invalid_type]), collapse = ", "),
            ". Treating as mismatch-only ('X').")
    hits$bulge_type[invalid_type] <- "X"
  }

  # ---- 8b. Deduplicate indexed-mode multi-representation hits ----
  # CRISPRitz -index mode writes all hit types (mismatch, DNA bulge, RNA bulge)
  # to .targets.txt without the merge step that normally produces .bestMerge.txt.
  # The same genomic locus may appear multiple times — once per hit type that
  # reaches it within the search parameters. Without deduplication, these
  # redundant representations inflate CFD sums and deflate specificity scores.
  #
  # Strategy: for each unique (guide, chr, pos, strand) locus, keep only the
  # representation with the lowest total edit distance. On ties, prefer
  # mismatch-only (bulge_type "X") over bulge hits, since CFD scoring is
  # most accurate for pure mismatch alignments.
  #
  # This is safe to run unconditionally — if the input has no duplicates

  # (e.g. from non-indexed .bestMerge.txt), the check is O(n) and exits
  # immediately.

  if (nrow(hits) > 0) {

    # Ensure total_edits is populated for sort ordering
    if (all(is.na(hits$total_edits))) {
      hits$total_edits <- hits$n_mismatches + hits$bulge_size
    } else {
      # Fill any sporadic NAs
      na_te <- is.na(hits$total_edits)
      if (any(na_te)) {
        hits$total_edits[na_te] <- hits$n_mismatches[na_te] + hits$bulge_size[na_te]
      }
    }

    # Build dedup key: one entry per (guide, locus)
    dedup_key <- paste(hits$guide_seq, hits$chr, hits$pos, hits$strand, sep = "\t")

    n_dupes <- sum(duplicated(dedup_key))

    if (n_dupes > 0) {
      n_before <- nrow(hits)

      # Sort to place the preferred representation first within each group:
      #   1. Lowest total_edits (fewest edits = most similar to guide = highest CFD)
      #   2. Lowest bulge_size (prefer pure mismatches over bulge alignments)
      #   3. bulge_type "X" first (explicit mismatch preference for CFD accuracy)
      bulge_type_rank <- ifelse(hits$bulge_type == "X", 0L,
                                ifelse(hits$bulge_type == "DNA", 1L,
                                       ifelse(hits$bulge_type == "RNA", 2L, 3L)))

      ord <- order(dedup_key, hits$total_edits, hits$bulge_size, bulge_type_rank,
                   na.last = TRUE)

      hits <- hits[ord, , drop = FALSE]
      dedup_key_sorted <- dedup_key[ord]

      keep <- !duplicated(dedup_key_sorted)
      hits <- hits[keep, , drop = FALSE]

      n_after <- nrow(hits)
      n_removed <- n_before - n_after
      pct_removed <- round(100 * n_removed / n_before, 1)

      # Diagnostic: report dedup summary (useful for validating indexed output)
      message("Deduplicated CRISPRitz indexed output: ",
              format(n_before, big.mark = ","), " -> ",
              format(n_after, big.mark = ","), " unique loci (",
              format(n_removed, big.mark = ","),
              " redundant representations removed, ", pct_removed, "%)")
    }
  }

  # ---- 9. Assemble output ----
  result <- data.frame(
    guide_seq       = hits$guide_seq,
    grna_indices    = I(hits$grna_indices),
    chr             = hits$chr,
    pos             = hits$pos,
    strand          = hits$strand,
    offtarget_seq   = hits$offtarget_seq_clean,
    pam_gen         = hits$pam_gen,
    n_mismatches    = hits$n_mismatches,
    bulge_type      = hits$bulge_type,
    bulge_size      = hits$bulge_size,
    total_edits     = hits$total_edits,
    stringsAsFactors = FALSE
  )

  # Retain gapped sequences for detailed reporting
  result$crRNA_aligned     <- hits$crRNA_raw
  result$offtarget_aligned <- hits$offtarget_seq

  return(result)
}
