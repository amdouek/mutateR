#' Design gRNA pairs for exon deletions
#'
#' Integrates exon structure, Cas effector targeting rules, exon phase compatibility,
#' protein domain annotations, and PTC predictions to propose candidate
#' gRNA pairs that could induce an in-frame (or allowable frame-shifting) exon deletion.
#'
#' @param gene_id Character. Gene identifier: either a gene symbol
#'        (via external_gene_name) or an Ensembl Gene ID (ENSG...).
#' @param species Character. Ensembl short code (e.g., "hsapiens").
#' @param id_type Character. One of "symbol" (default) or "ensembl_gene_id".
#' @param nuclease Character. Must be specified explicitly. 
#'        Currently supported: "Cas9", "Cas12a".
#' @param filter_domains Character vector. Optional list of InterPro IDs or descriptions
#'                       to focus designs onto (e.g., c("IPR008967", "DNA-binding")).
#' @param transcript_id Character. Optional transcript ID to override canonical 
#'                      transcript choice. Warns if not canonical.
#' @param genome BSgenome. Required BSgenome object for target species.
#' @param return_all Logical (default FALSE). If TRUE, will retrieve ALL exon pairs 
#'                   (both phase compatible and incompatible). This is potentially huge, so interactive confirmation is required.
#'
#' @return A data.frame of candidate gRNA pairs with columns:
#'   - transcript_id
#'   - exon_5p / exon_3p
#'   - gRNA_5p / pam_5p
#'   - gRNA_3p / pam_3p
#'   - compatible (exon phase check TRUE/FALSE)
#'   - frameshift (TRUE if not divisible by 3)
#'   - ptc_flag (TRUE if frameshift causes PTC in non-terminal exon)
#'   - terminal_exon_case (TRUE if frameshift causes PTC in terminal exon)
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' results <- design_grna_pairs("TP53", "hsapiens",
#'                              nuclease="Cas9",
#'                              genome=BSgenome.Hsapiens.UCSC.hg38)
#' head(results)
#' }
#' @export

design_grna_pairs <- function(gene_id, species,
                              id_type = c("symbol", "ensembl_gene_id"),
                              nuclease, 
                              filter_domains = NULL, 
                              transcript_id = NULL,
                              genome,
                              return_all = FALSE) {
  
  id_type <- match.arg(id_type)
  
  # ----- Supported nucleases -----
  valid_nucleases <- c("Cas9", "Cas12a")
  if (missing(nuclease)) {
    stop("You must explicitly specify a nuclease (e.g. nuclease='Cas9' or nuclease='Cas12a').")
  }
  if (!nuclease %in% valid_nucleases) {
    stop("Unsupported nuclease '", nuclease, "'. Allowed: ", paste(valid_nucleases, collapse=", "))
  }
  
  # ------ Step 1: Gene/transcript info ------
  tx_info <- get_gene_info(gene_id, species, id_type = id_type)
  if (is.null(tx_info$all) || nrow(tx_info$all) == 0) {
    stop("No transcript found for gene: ", gene_id)
  }
  
  if (is.null(transcript_id)) {
    if (!is.null(tx_info$canonical)) {
      tx_id <- tx_info$canonical$ensembl_transcript_id[1]
    } else {
      warning("No canonical transcript annotated for ", gene_id, 
              "; defaulting to first transcript returned.")
      tx_id <- tx_info$all$ensembl_transcript_id[1]
    }
  } else {
    tx_id <- transcript_id
    if (tx_id %in% tx_info$all$ensembl_transcript_id) {
      if (!tx_id %in% tx_info$canonical$ensembl_transcript_id) {
        warning("Selected transcript ", tx_id, 
                " is NOT the canonical transcript for ", gene_id)
      }
    } else {
      warning("Transcript ", tx_id, " not found for ", gene_id, " (species: ", species, ").")
    }
  }
  message("Using transcript: ", tx_id)
  
  # ------ Step 2: Exon structures ------
  exons <- get_exon_structures(tx_id, species, output = "GRanges")
  exon_meta <- as.data.frame(mcols(exons))
  exon_meta$rank <- seq_len(nrow(exon_meta))
  
  ## ------ Handling for ≤2 total or coding exons ------
  n_total_exons  <- nrow(exon_meta)
  n_coding_exons <- sum(!is.na(exon_meta$cds_start) & !is.na(exon_meta$cds_end))
  
  if (n_total_exons <= 2 || n_coding_exons <= 2) {
    warning("Transcript ", tx_id, " has ", n_total_exons, " total exon(s) (",
            n_coding_exons, " coding). Performing intragenic large‑deletion design.")
    
    ### ----- Step 2A: Locate gRNA sites (Cas9 only for now) -----
    if (nuclease != "Cas9")
      warning("Edge‑case mode currently implemented for Cas9 only.")
    
    all_sites <- find_cas9_sites(exons, genome, pam = "NGG", protospacer_length = 20)
    if (is.null(all_sites) || length(all_sites) == 0) {
      warning("No Cas9 sites found in the provided exons.")
      return(list(transcript_id = tx_id,
                  exons = exons,
                  scored_grnas = NULL,
                  pairs = data.frame(),
                  plot = plot_grna_design(exons, NULL, transcript_id = tx_id)))
    }
    
    ## ----- Step 2B: Score and filter guides -----
    scored_sites <- score_grnas(all_sites, method = "ruleset1")
    on_target <- as.numeric(mcols(scored_sites)$ontarget_score)
    min_cutoff <- 0.4
    keep <- which(!is.na(on_target) & on_target >= min_cutoff)
    filtered_sites <- scored_sites[keep]
    message(length(filtered_sites), " gRNAs above on‑target score cutoff (", min_cutoff, ").")
    
    if (length(filtered_sites) < 2) {
      warning("Fewer than two high‑scoring guides; deletion design not possible.")
      return(list(transcript_id = tx_id,
                  exons = exons,
                  scored_grnas = filtered_sites,
                  pairs = data.frame(),
                  plot = plot_grna_design(exons, NULL, transcript_id = tx_id)))
    }
    
    ## ----- Step 2C: create all possible intragenic pairs -----
    site_df <- as.data.frame(filtered_sites)
    site_df$pos <- GenomicRanges::start(filtered_sites)
    comb_idx <- utils::combn(seq_len(nrow(site_df)), 2)
    pair_df <- data.frame(
      idx5  = comb_idx[1,],
      idx3  = comb_idx[2,],
      exon_5p = site_df$exon_rank[comb_idx[1,]],
      exon_3p = site_df$exon_rank[comb_idx[2,]],
      start_5p = site_df$pos[comb_idx[1,]],
      start_3p = site_df$pos[comb_idx[2,]],
      seq_5p = site_df$protospacer_sequence[comb_idx[1,]],
      seq_3p = site_df$protospacer_sequence[comb_idx[2,]],
      score_5p = on_target[comb_idx[1,]],
      score_3p = on_target[comb_idx[2,]],
      del_size = abs(site_df$pos[comb_idx[2,]] - site_df$pos[comb_idx[1,]])
    )
    
    pair_df <- pair_df[order(-pair_df$del_size), ]
    top_n <- min(50, nrow(pair_df))
    pair_df <- utils::head(pair_df, top_n)
    
    score_cutoff <- 0.6
    pair_df$recommended <- with(pair_df,
                                score_5p >= score_cutoff & score_3p >= score_cutoff)
    
    message("Generated ", nrow(pair_df),
            " intragenic deletion pairs; ",
            sum(pair_df$recommended), " meet score cutoff.")
    
    return(list(
      transcript_id = tx_id,
      exons = exons,
      scored_grnas = filtered_sites,
      pairs = pair_df,
      plot = plot_grna_design(exons, NULL, transcript_id = tx_id)
    ))
  } else {
    # guard to stop function after early return in edge‑case mode
    # (needed in case of single‑exon transcript so code below never runs)
    if (nrow(as.data.frame(mcols(exons))) <= 2) return(invisible(NULL))
  }
  
  
  # ----- Step 3: Find gRNA sites and retain genomic coordinates -----
  if (nuclease == "Cas9") {
    grna_sites <- find_cas9_sites(exons, genome, pam = "NGG", protospacer_length = 20)
  } else {
    grna_sites <- find_cas12a_sites(exons, genome, pam = "TTTV", protospacer_length = 23)
  }
  if (is.null(grna_sites)) {
    warning("No gRNA sites found for ", gene_id, " transcript ", tx_id)
    return(NULL)
  }
  
  ## ----- Attach exon rank if needed -----
  if (!"exon_rank" %in% names(mcols(grna_sites))) {
    overlaps <- GenomicRanges::findOverlaps(grna_sites, exons)
    grna_sites$exon_rank <- NA_integer_
    grna_sites$exon_rank[queryHits(overlaps)] <- exon_meta$rank[subjectHits(overlaps)]
  }
  
  ## ----- Extract as data.frame for easy merging -----
  grna_meta <- as.data.frame(grna_sites)
  grna_meta$seqnames <- as.character(GenomicRanges::seqnames(grna_sites))
  grna_meta$start <- GenomicRanges::start(grna_sites)
  grna_meta$end <- GenomicRanges::end(grna_sites)
  grna_meta$strand <- as.character(GenomicRanges::strand(grna_sites))
  
  # ----- Step 4: Phase compatible exon pairs -----
  comp_pairs <- check_exon_phase(exon_meta)
  if (!return_all) {
    comp_pairs <- subset(comp_pairs, compatible == TRUE)
  } else {
    if (interactive()) {
      ans <- tryCatch({
        readline("You requested gRNAs from ALL exon pairs (including phase-incompatible ones). \n This may generate an enormous number of gRNA pairs. \n Are you sure you want to proceed? [y/N]: ")
      }, interrupt = function(e) {
        cat("\nInterrupted. Aborting.\n"); return("n")
      })
      if (tolower(ans) != "y") stop("Aborted by user.")
      else message("Proceeding with ALL exon pairs (warning: large output possible).")
    } else {
      stop("return_all=TRUE disallowed in non-interactive mode.")
    }
  }
  
  # ----- Step 5: Map InterPro domains -----
  domain_map <- map_protein_domains(tx_id, species)
  
  # ----- Step 6: Filter by domain if requested -----
  if (!is.null(filter_domains) && !is.null(domain_map)) {
    domain_map <- subset(domain_map, domain_id %in% filter_domains |
                           domain_desc %in% filter_domains)
    if (nrow(domain_map) == 0) warning("No InterPro domains found matching filter criteria.")
    exons_to_consider <- unique(domain_map$exon_rank)
    comp_pairs <- comp_pairs[comp_pairs$exon_5p %in% exons_to_consider |
                               comp_pairs$exon_3p %in% exons_to_consider, ]
  }
  
  # ----- Step 7: Build gRNA pairs -----
  results <- list()
  for (i in seq_len(nrow(comp_pairs))) {
    e5 <- comp_pairs$exon_5p[i]
    e3 <- comp_pairs$exon_3p[i]
    
    g5 <- subset(grna_meta, exon_rank == e5)
    g3 <- subset(grna_meta, exon_rank == e3)
    
    if (nrow(g5) > 0 && nrow(g3) > 0) {
      combos <- merge(g5, g3, by = NULL, suffixes = c("_5p", "_3p"))
      combos$exon_5p <- e5
      combos$exon_3p <- e3
      combos$compatible <- comp_pairs$compatible[i]
      combos$transcript_id <- tx_id
      
      ## ----- Frameshift/PTC check -----
      fs_check <- check_frameshift_ptc(exon_meta, e5, e3)
      combos$frameshift         <- fs_check$frameshift
      combos$ptc_flag           <- fs_check$ptc_flag
      combos$terminal_exon_case <- fs_check$terminal_exon_case
      
      results[[length(results) + 1]] <- combos
    }
  }
  
  if (length(results) == 0) {
    warning("No valid gRNA pairs identified.")
    return(NULL)
  }
  
  out <- do.call(rbind, results)
  
  # ----- Step 8: Select and order relevant columns cleanly -----
  out <- out[, c(
    "transcript_id",
    "exon_5p", "exon_3p", "compatible",
    "protospacer_sequence_5p", "pam_sequence_5p",
    "protospacer_sequence_3p", "pam_sequence_3p",
    "seqnames_5p", "start_5p", "end_5p", "strand_5p",
    "seqnames_3p", "start_3p", "end_3p", "strand_3p",
    "frameshift", "ptc_flag", "terminal_exon_case"
  )]
  
  return(out)
}