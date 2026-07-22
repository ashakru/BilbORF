#' Annotate ORFs with Translation Status Across Transcript Isoforms
#'
#' Maps Open Reading Frames (ORFs) to overlapping transcript isoforms and
#' determines their translation competency by analyzing start codons, stop codons,
#' and internal stops. Each ORF-transcript pair is evaluated independently.
#' Every input ORF appears in the output at least once.
#'
#' @param annotations List containing genomic annotations, as returned by
#'   \code{\link{prepare_annotations_fromGTF}}. Must contain a \code{transcripts}
#'   element with transcript models as a GRangesList.
#' @param orfs GRanges or GRangesList with genomic coordinates of ORF start and
#'   stop codons. For GRanges input, ranges will be grouped by name. Each ORF
#'   can span multiple exons (for spliced ORFs).
#' @param BSgenome A BSgenome object matching the genome build used for annotations.
#'   Required for extracting ORF nucleotide sequences.
#' @param transcript_meta A data.frame with transcript metadata. Must contain
#'   a \code{transcript_id} column matching \code{names(annotations$transcripts)}.
#'   Commonly includes: gene_name, gene_id, transcript_type. This information
#'   is joined to the output table.
#' @param orfs_meta Optional data.frame with ORF metadata. Must contain an
#'   \code{ORF_id} column matching ORF names. Additional columns are joined
#'   to the output.
#' @param start_codons Character vector of valid start codons. Default is
#'   \code{c("ATG","TTG","CTG","GTG")} covering canonical and alternative starts.
#' @param stop_codons Character vector of valid stop codons. Default is
#'   \code{c("TAG", "TAA", "TGA")} covering all standard stop codons.
#' @param stop_codon_convention How the input 3' boundary represents the stop
#'   codon. `"included"` validates the terminal triplet, `"excluded"` validates
#'   and adds the following three transcript-oriented nucleotides, and
#'   `"auto"` accepts a convention only when exactly one candidate is a valid
#'   stop. The clean parsed-caller cache should use `"excluded"` explicitly.
#' @param check_stop_codon Deprecated compatibility argument. The 3 nt
#'   immediately 3' of each input ORF are now always checked because they are
#'   required to resolve and validate the stop-codon convention. Splice
#'   junctions are handled correctly: if the ORF ends in one exon the search
#'   continues into the next exon. The output includes:
#'   \code{downstream_codon} (the 3 nt, or \code{NA} if the transcript ends
#'   within 3 nt), \code{downstream_is_stop} (logical), and
#'   \code{stop_codon_end} (strand-aware genomic position of the last nt of
#'   the downstream triplet: \code{max(end)} for + strand, \code{min(start)}
#'   for - strand). The argument currently has no effect. Default \code{FALSE}.
#' @param cds_gr Optional \code{GRangesList} of annotated CDS exons, keyed by
#'   \code{transcript_id} (the same key space as
#'   \code{names(annotations\$transcripts)}). Typically obtained via
#'   \code{cdsBy(txdb, by = "tx", use.names = TRUE)} or by subsetting a
#'   pre-built list to the relevant transcripts. When supplied, each ORF-
#'   transcript pair is tested for exon-level overlap with the CDS of that
#'   specific transcript, and the reading frame is determined in transcript
#'   coordinates (so splice junctions between the ORF and the CDS are handled
#'   correctly). Frame is \code{(tx_coord_orf_5prime - tx_coord_cds_5prime)
#'   \%\% 3 == 0}. Adds \code{overlaps_cds} (logical, \code{FALSE} if no
#'   overlap, \code{NA} for intronic / no-transcript rows) and
#'   \code{cds_frame} (\code{"in_frame"} / \code{"out_of_frame"} /
#'   \code{NA}). It also classifies every successfully normalised ORF-isoform
#'   pair into \code{reference_orf_type} and \code{reference_orf_class} using
#'   the already-computed transcript coordinates; no second ORF/transcript
#'   intersection is performed. Default \code{NULL}.
#' @param chunk_size Integer. When set, ORFs are split into batches of this
#'   size and processed sequentially, with \code{gc()} called between batches
#'   to release memory. This keeps peak memory proportional to
#'   \code{chunk_size} rather than the full dataset. \code{unique_aa_id} is
#'   recomputed globally after all chunks are combined. A value of 500-2000 is
#'   usually appropriate. Note: true parallelism is not supported because
#'   \code{BSgenome} objects contain C-level pointers that do not serialise
#'   safely across worker processes. Default \code{NULL} (no chunking).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{table}{A data.frame (tibble) with one row per ORF-transcript pair.
#'       Every input ORF appears at least once. Key columns include:
#'       \itemize{
#'         \item \code{ORF_isoform_id}: Unique identifier for ORF-transcript pair
#'         \item \code{ORF_id}: Original ORF identifier
#'         \item \code{transcript_id}: Overlapping transcript ID (NA for
#'           \code{orf_status = "no_transcript"})
#'         \item \code{gene_id}, \code{gene_name}: Gene annotations
#'         \item \code{seq_nt}: Nucleotide sequence of the ORF in this isoform
#'           (NA for \code{"intronic"} and \code{"no_transcript"})
#'         \item \code{seq_aa}: Translated amino acid sequence
#'         \item \code{len_nt}, \code{len_aa}: Sequence lengths
#'         \item \code{start_codon}, \code{stop_codon}: First and last codons
#'         \item \code{orf_status}: Translation status (see Details)
#'         \item \code{complete_codons}: FALSE if sequence length is not
#'           divisible by 3; NA for rows without sequence
#'         \item \code{unique_tx_iso}: Groups isoforms producing identical
#'           proteins (NA for non-translatable statuses)
#'       }}
#'     \item{ranges}{A stop-inclusive `GRangesList` of successfully projected
#'       ORF-transcript pairs. Names match `ORF_isoform_id` in `table`. The table
#'       reports both stop-exclusive and stop-inclusive genomic bounds.}
#'   }
#'
#' @details
#' \strong{ORF Status Categories:}
#' \describe{
#'   \item{translatable}{Valid start and stop codons present, no internal stops.
#'     The ORF can produce a full-length protein in this isoform.}
#'   \item{internal_stop}{Contains one or more stop codons before the final
#'     position. Translation would terminate prematurely.}
#'   \item{no_stop}{Missing a valid stop codon at the 3' end. ORF extends to
#'     transcript terminus without proper termination.}
#'   \item{no_start}{Missing a valid start codon at the 5' end but has a stop.
#'     Cannot initiate translation properly.}
#'   \item{no_stop_no_start}{Missing both valid start and stop codons.}
#'   \item{intronic}{The ORF overlaps a transcript's genomic span but all its
#'     exons fall within introns of that transcript. Sequence columns are NA.}
#'   \item{no_transcript}{The ORF does not overlap any annotated transcript.
#'     Transcript and sequence columns are NA.}
#' }
#'
#' Out-of-frame ORF-isoform pairs (where sequence length is not divisible by 3)
#' are flagged via \code{complete_codons = FALSE} and a warning is emitted, but
#' they are still classified with the standard status categories.
#'
#' \strong{Algorithm Overview:}
#' \enumerate{
#'   \item Find overlaps between ORFs and transcripts (GRangesList level);
#'     ORFs with no overlap are recorded as \code{"no_transcript"}
#'   \item Vectorised exon-level intersection: unlist both GRLs, single
#'     \code{findOverlaps} + \code{pintersect} + \code{splitAsList};
#'     entirely intronic pairs are recorded as \code{"intronic"}
#'   \item Extract and translate sequences for the remaining pairs
#'   \item Classify translation status
#'   \item Group isoforms producing identical proteins
#'   \item Append stub rows for intronic and no-transcript ORFs
#' }
#'
#' @note Set \code{stop_codon_convention} to match whether the input coordinates
#'   contain the terminal stop codon. Irrespective of the input convention,
#'   \code{ranges} and the primary \code{start}/\code{end}/sequence columns are
#'   stop-inclusive. For multi-exon ORFs, provide all exonic segments.
#'
#' @importFrom ORFik groupGRangesBy
#' @import tibble
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom IRanges relist PartitioningByEnd
#' @importFrom dplyr mutate case_when left_join group_by ungroup relocate bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(rtracklayer)
#'
#' # Prepare annotations
#' gtf <- "path/to/annotation.gtf"
#' annotations <- prepare_annotations_fromGTF(gtf, BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Load ORFs from BED file
#' orfs <- import("path/to/orfs.bed", format = "BED")
#' names(orfs) <- orfs$name
#'
#' # Prepare transcript metadata
#' transcripts_meta <- import(gtf, format = "GTF") %>%
#'   as.data.frame() %>%
#'   filter(type == "transcript") %>%
#'   select(gene_name, gene_id, transcript_id, transcript_type) %>%
#'   distinct()
#'
#' # Annotate ORFs
#' result <- annotate_orf_isoforms(
#'   annotations,
#'   orfs,
#'   BSgenome.Hsapiens.UCSC.hg38,
#'   transcripts_meta
#' )
#'
#' # Explore results
#' table(result$table$orf_status)
#' head(result$table)
#' }
#'
#' # Using package example data
#' gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf",
#'                    package = "BilbORF")
#' bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
#'
#' if (gtf != "" && bed != "") {
#'   library(BSgenome.Hsapiens.UCSC.hg38)
#'   library(rtracklayer)
#'
#'   BSgenome <- BSgenome.Hsapiens.UCSC.hg38
#'   annotations <- prepare_annotations_fromGTF(gtf, BSgenome)
#'   orfs <- import(bed, format = "BED")
#'   names(orfs) <- orfs$name
#'
#'   transcripts_meta <- import(gtf, format = "GTF") %>%
#'     as.data.frame() %>%
#'     dplyr::filter(type == "transcript") %>%
#'     dplyr::select(gene_name, gene_id, transcript_id, transcript_type) %>%
#'     dplyr::distinct()
#'
#'   annotated_orfs <- annotate_orf_isoforms(
#'     annotations, orfs[1:5], BSgenome, transcripts_meta
#'   )
#'
#'   print(table(annotated_orfs$table$orf_status))
#' }
#'
#' @seealso
#' \code{\link{prepare_annotations_fromGTF}} for creating annotation input
#' \code{\link{diff_orf_usage}} for differential usage analysis
annotate_orf_isoforms <- function(annotations, orfs, BSgenome,
                                  transcript_meta, orfs_meta = NULL,
                                  start_codons     = c("ATG","TTG","CTG","GTG"),
                                  stop_codons      = c("TAG", "TAA", "TGA"),
                                  stop_codon_convention = c(
                                    "auto", "included", "excluded"
                                  ),
                                  check_stop_codon = FALSE,
                                  cds_gr           = NULL,
                                  chunk_size       = NULL) {

  stop_codon_convention <- match.arg(stop_codon_convention)

  # --- Chunked processing ---
  # When chunk_size is set and there are more ORFs than chunk_size, split the
  # input into batches, process each independently, then combine.  This keeps
  # peak memory proportional to chunk_size rather than the full dataset size.
  if (!is.null(chunk_size)) {
    if (!is(orfs, "GRangesList")) orfs <- ORFik::groupGRangesBy(orfs)
    n <- length(orfs)
    if (n > chunk_size) {
      chunks   <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
      n_chunks <- length(chunks)
      cli::cli_inform(c(
        "i" = "Processing {n} ORFs in {n_chunks} chunk{?s} of up to {chunk_size}."
      ))
      results <- vector("list", n_chunks)
      for (k in seq_along(chunks)) {
        cli::cli_inform("  Chunk {k}/{n_chunks} ({length(chunks[[k]])} ORFs)")
        results[[k]] <- annotate_orf_isoforms(
          annotations      = annotations,
          orfs             = orfs[chunks[[k]]],
          BSgenome         = BSgenome,
          transcript_meta  = transcript_meta,
          orfs_meta        = orfs_meta,
          start_codons     = start_codons,
          stop_codons      = stop_codons,
          stop_codon_convention = stop_codon_convention,
          check_stop_codon = check_stop_codon,
          cds_gr           = cds_gr,
          chunk_size       = NULL    # prevent recursive chunking
        )
        gc(verbose = FALSE)          # release intermediate objects before next chunk
      }
      # Combine tables and recompute unique_aa_id globally
      combined_table <- dplyr::bind_rows(lapply(results, `[[`, "table"))
      unique_seqs <- unique(combined_table$seq_aa[!is.na(combined_table$seq_aa)])
      combined_table$unique_aa_id <- dplyr::if_else(
        is.na(combined_table$seq_aa),
        NA_character_,
        paste0("aa_iso_", match(combined_table$seq_aa, unique_seqs))
      )
      combined_ranges <- do.call(c, lapply(results, `[[`, "ranges"))
      return(list(table = combined_table, ranges = combined_ranges))
    }
  }

  tx <- annotations$transcripts

  # --- Input validation ---
  if (is.null(orfs)) {
    cli::cli_abort("{.arg orfs} must not be NULL")
  }

  if (!is(orfs, "GRangesList")) {
    orfs <- ORFik::groupGRangesBy(orfs)
  }

  if (length(orfs) == 0) {
    cli::cli_abort("{.arg orfs} must not be empty")
  }

  # Helper: min/max genomic bounds per element of a GRangesList.
  # range(grl) applies range() elementwise (C-level), giving one bounding
  # GRanges per list element; unlist flattens to a flat GRanges.
  grl_bounds <- function(grl) {
    rng <- unlist(range(grl), use.names = FALSE)
    list(
      seqnames = as.character(seqnames(rng)),
      start    = as.integer(start(rng)),
      end      = as.integer(end(rng))
    )
  }

  # --- Find ORF-transcript overlaps (GRangesList level) ---
  orf_tx_overlap <- findOverlaps(orfs, tx)

  # Track ORFs that don't overlap any transcript
  unmapped_orf_idx <- setdiff(seq_along(orfs), unique(queryHits(orf_tx_overlap)))

  n_pairs   <- length(orf_tx_overlap)
  orf_by_tx <- orfs[queryHits(orf_tx_overlap)]
  tx_by_orf <- tx[subjectHits(orf_tx_overlap)]

  # --- Vectorised pairwise exon intersection ---
  # Unlist both GRLs to flat GRanges, tracking pair membership.
  # A single findOverlaps + pintersect replaces per-pair mapply calls.
  orf_exons <- unlist(orf_by_tx, use.names = FALSE)
  tx_exons  <- unlist(tx_by_orf, use.names = FALSE)
  mcols(orf_exons) <- NULL
  mcols(tx_exons)  <- NULL

  orf_pair <- rep(seq_len(n_pairs), elementNROWS(orf_by_tx))
  tx_pair  <- rep(seq_len(n_pairs), elementNROWS(tx_by_orf))

  # Find exon-level overlaps; retain only hits within the same ORF-tx pair
  exon_hits <- findOverlaps(orf_exons, tx_exons)
  same_pair <- orf_pair[queryHits(exon_hits)] == tx_pair[subjectHits(exon_hits)]
  exon_hits <- exon_hits[same_pair]

  isect      <- pintersect(orf_exons[queryHits(exon_hits)],
                           tx_exons[subjectHits(exon_hits)])
  mcols(isect)  <- NULL  # strip 'hit' metadata added by pintersect
  isect_pair    <- orf_pair[queryHits(exon_hits)]

  # Rebuild GRangesList indexed 1..n_pairs.
  # relist + PartitioningByEnd is robust and always returns exactly n_pairs
  # elements (empty slots for intronic pairs).
  # Requires isect_pair to be non-decreasing — guaranteed because findOverlaps
  # returns hits sorted by queryHits, and orf_pair is non-decreasing.
  counts_per_pair <- tabulate(isect_pair, nbins = n_pairs)
  orf_in_tx_all   <- relist(isect, PartitioningByEnd(cumsum(counts_per_pair)))

  # Separate intronic pairs from pairs with exonic overlap
  keep         <- which(lengths(orf_in_tx_all) > 0)
  intronic_idx <- which(lengths(orf_in_tx_all) == 0)

  intronic_orf_ids  <- names(orf_by_tx)[intronic_idx]
  intronic_tx_ids   <- names(tx_by_orf)[intronic_idx]
  intronic_orfs_grl <- orf_by_tx[intronic_idx]

  orf_in_tx <- orf_in_tx_all[keep]
  orf_by_tx <- orf_by_tx[keep]
  tx_by_orf <- tx_by_orf[keep]

  # --- Rank exons 5' to 3' ---
  orf_in_tx <- rank_exons(orf_in_tx)

  # --- Set ORF isoform names ---
  ORF_isoform_id <- paste(names(orf_by_tx), paste0("iso-", names(tx_by_orf)), sep = "_")
  names(orf_in_tx) <- ORF_isoform_id

  # --- Extract and translate sequences ---
  seq_nt     <- extractTranscriptSeqs(BSgenome, orf_in_tx)
  seq_aa     <- suppressWarnings(Biostrings::translate(seq_nt))
  seq_nt_chr <- as.character(seq_nt)
  seq_aa_chr <- as.character(seq_aa)

  # --- Vectorised genomic coordinates for mapped pairs ---
  bc         <- grl_bounds(orf_in_tx)
  genomic_sn <- bc$seqnames
  genomic_st <- bc$start
  genomic_en <- bc$end

  # Strand for each ORF-tx pair: unlist once and index into the first exon of
  # each element (all exons in a pair share the same strand).
  strands_vec <- suppressWarnings({
    ul_s      <- as.character(strand(unlist(orf_in_tx, use.names = FALSE)))
    first_idx <- cumsum(c(1L, elementNROWS(orf_in_tx)[-length(orf_in_tx)]))
    ul_s[first_idx]
  })

  # --- Transcript-coordinate stop candidates ---
  downstream_codon_vec   <- rep(NA_character_, length(orf_in_tx))
  downstream_is_stop_vec <- rep(NA,             length(orf_in_tx))
  stop_codon_end_vec     <- rep(NA_integer_,    length(orf_in_tx))

  {
    # --- Vectorised downstream stop codon check ---
    # Key insight: extract full transcript sequence once per unique transcript,
    # then downstream codon = substr(tx_seq, orf_end_in_tx + 1, orf_end_in_tx + 3).
    # Cost: one extractTranscriptSeqs call + O(n_unique_tx) coordinate work.

    # 1. Deduplicated transcript GRangesList (one entry per unique tx ID)
    unique_tx_names <- unique(names(tx_by_orf))
    first_occ       <- match(unique_tx_names, names(tx_by_orf))
    tx_grl_unique   <- tx_by_orf[first_occ]
    names(tx_grl_unique) <- unique_tx_names

    # 2. Full transcript sequences — single BSgenome query
    tx_seqs_full <- suppressWarnings(
      as.character(extractTranscriptSeqs(BSgenome, tx_grl_unique))
    )
    names(tx_seqs_full) <- unique_tx_names

    # 3. Vectorised ORF 5' position in transcript coordinates.
    #    Use split() for O(n) grouping instead of which()==tid per transcript.
    orf_5p_g_all  <- ifelse(strands_vec == "+", genomic_st, genomic_en)
    orf_5p_tx_all <- rep(NA_integer_, length(orf_in_tx))
    tx_groups     <- split(seq_along(orf_in_tx), names(tx_by_orf))

    for (tid in names(tx_groups)) {
      idx   <- tx_groups[[tid]]
      strnd <- strands_vec[idx[1L]]
      tx_ex <- tx_grl_unique[[tid]]
      mcols(tx_ex) <- NULL
      orf_5p_tx_all[idx] <- .genomic_5prime_to_tx_coord_vec(
        orf_5p_g_all[idx], tx_ex, strnd
      )
    }

    # 4. substr-based downstream codon (1-based R indexing)
    #    orf_5p_tx_all is 0-based; ORF occupies positions [orf_5p_tx, orf_5p_tx+len)
    #    downstream starts at position orf_5p_tx + len (0-based) = +1 in 1-based
    len_nt_check  <- nchar(seq_nt_chr)
    ds_start_1b   <- orf_5p_tx_all + len_nt_check + 1L   # 1-based start
    ds_end_1b     <- ds_start_1b + 2L                    # 1-based end (3 nt)

    tx_seq_per_orf <- tx_seqs_full[names(tx_by_orf)]
    tx_len_per_orf <- nchar(tx_seq_per_orf)
    valid_ds       <- !is.na(orf_5p_tx_all) & ds_end_1b <= tx_len_per_orf

    if (any(valid_ds)) {
      downstream_codon_vec[valid_ds] <- substr(
        tx_seq_per_orf[valid_ds], ds_start_1b[valid_ds], ds_end_1b[valid_ds]
      )
      downstream_is_stop_vec[valid_ds] <-
        downstream_codon_vec[valid_ds] %in% stop_codons
    }

    # 5. stop_codon_end: genomic coordinate of the 3'-most base of the
    #    downstream codon.  Convert 0-based tx coord (ds_end_1b - 1) back to
    #    genomic using sorted exon structure — no BSgenome access needed.
    ds_end_0b <- ds_end_1b - 1L   # 0-based position of last ds nt in tx

    for (tid in names(tx_groups)) {
      idx        <- tx_groups[[tid]]
      valid_here <- idx[valid_ds[idx]]
      if (length(valid_here) == 0L) next

      strnd <- strands_vec[valid_here[1L]]
      tx_ex <- tx_grl_unique[[tid]]
      mcols(tx_ex) <- NULL
      tx_s  <- if (strnd == "+") sort(tx_ex) else rev(sort(tx_ex))
      ex_w  <- as.integer(width(tx_s))
      cum_b <- c(0L, cumsum(ex_w)[-length(ex_w)])

      pos <- ds_end_0b[valid_here]
      k   <- findInterval(pos, cum_b)
      ok  <- k >= 1L & k <= length(tx_s)

      if (strnd == "+") {
        stop_codon_end_vec[valid_here[ok]] <-
          as.integer(start(tx_s))[k[ok]] + (pos[ok] - cum_b[k[ok]])
      } else {
        stop_codon_end_vec[valid_here[ok]] <-
          as.integer(end(tx_s))[k[ok]] - (pos[ok] - cum_b[k[ok]])
      }
    }
  }

  # --- Resolve convention and construct stop-inclusive output chains ---
  # `orf_5p_tx_all` and the transcript sequences above are reused; no second
  # ORF/transcript exon intersection is performed.
  input_seq_nt_chr <- seq_nt_chr
  input_len_nt <- nchar(input_seq_nt_chr)
  terminal_codon_vec <- ifelse(
    input_len_nt >= 3L,
    substr(input_seq_nt_chr, input_len_nt - 2L, input_len_nt),
    NA_character_
  )
  terminal_is_stop_vec <- terminal_codon_vec %in% stop_codons

  resolved_stop_convention <- rep(NA_character_, length(orf_in_tx))
  stop_adjustment_nt <- rep(NA_integer_, length(orf_in_tx))
  stop_normalization_status <- rep(
    "stop convention unresolved", length(orf_in_tx)
  )

  if (stop_codon_convention == "included") {
    resolved_stop_convention[] <- "included"
    stop_adjustment_nt[] <- 0L
    stop_normalization_status[] <- "input coordinates include stop codon"
  } else if (stop_codon_convention == "excluded") {
    resolved_stop_convention[] <- "excluded"
    stop_adjustment_nt[] <- 3L
    stop_normalization_status[] <-
      "input coordinates exclude stop codon; added 3 transcript nt"
  } else {
    included_only <- terminal_is_stop_vec & !(downstream_is_stop_vec %in% TRUE)
    excluded_only <- downstream_is_stop_vec %in% TRUE & !terminal_is_stop_vec
    both <- terminal_is_stop_vec & downstream_is_stop_vec %in% TRUE
    neither <- !terminal_is_stop_vec & !(downstream_is_stop_vec %in% TRUE)

    resolved_stop_convention[included_only] <- "included"
    stop_adjustment_nt[included_only] <- 0L
    stop_normalization_status[included_only] <-
      "auto: input coordinates include stop codon"
    resolved_stop_convention[excluded_only] <- "excluded"
    stop_adjustment_nt[excluded_only] <- 3L
    stop_normalization_status[excluded_only] <-
      "auto: input coordinates exclude stop codon; added 3 transcript nt"
    stop_normalization_status[both] <-
      "auto: terminal and downstream triplets are both stop codons"
    stop_normalization_status[neither] <-
      "auto: neither terminal nor downstream triplet is a stop codon"
  }

  input_3p_tx_0b <- orf_5p_tx_all + input_len_nt - 1L
  inclusive_3p_tx_0b <- input_3p_tx_0b + stop_adjustment_nt
  exclusive_3p_tx_0b <- inclusive_3p_tx_0b - 3L
  canonical_width_nt <- inclusive_3p_tx_0b - orf_5p_tx_all + 1L
  coordinate_projection_valid <-
    !is.na(orf_5p_tx_all) & !is.na(inclusive_3p_tx_0b) &
    exclusive_3p_tx_0b >= orf_5p_tx_all &
    inclusive_3p_tx_0b < tx_len_per_orf

  canonical_stop_codon <- ifelse(
    resolved_stop_convention == "included",
    terminal_codon_vec,
    ifelse(
      resolved_stop_convention == "excluded",
      downstream_codon_vec,
      NA_character_
    )
  )
  canonical_stop_valid <-
    coordinate_projection_valid & canonical_stop_codon %in% stop_codons
  canonical_complete_codons <-
    coordinate_projection_valid & canonical_width_nt %% 3L == 0L
  stop_normalization_valid <-
    canonical_stop_valid & canonical_complete_codons

  stop_normalization_status[
    !is.na(stop_adjustment_nt) & !coordinate_projection_valid
  ] <- "normalised stop lies outside the compatible transcript"
  stop_normalization_status[
    coordinate_projection_valid & !canonical_stop_valid
  ] <- "normalised terminal triplet is not a permitted stop codon"
  stop_normalization_status[
    canonical_stop_valid & !canonical_complete_codons
  ] <- "normalised ORF length is not divisible by three"

  projected_idx <- which(coordinate_projection_valid)
  ranges_stop_inclusive <- GenomicRanges::GRangesList(
    lapply(projected_idx, function(i) {
      .slice_transcript_interval(
        tx_by_orf[[i]],
        orf_5p_tx_all[i] + 1L,
        inclusive_3p_tx_0b[i] + 1L
      )
    })
  )
  names(ranges_stop_inclusive) <- ORF_isoform_id[projected_idx]

  ranges_stop_exclusive <- GenomicRanges::GRangesList(
    lapply(projected_idx, function(i) {
      .slice_transcript_interval(
        tx_by_orf[[i]],
        orf_5p_tx_all[i] + 1L,
        exclusive_3p_tx_0b[i] + 1L
      )
    })
  )
  names(ranges_stop_exclusive) <- ORF_isoform_id[projected_idx]

  pair_range_bounds <- function(paths, pair_index, n_pairs) {
    result <- list(
      start = rep(NA_integer_, n_pairs),
      end = rep(NA_integer_, n_pairs),
      orf_3p = rep(NA_integer_, n_pairs)
    )
    if (!length(paths)) return(result)
    bounds <- unlist(range(paths), use.names = FALSE)
    path_strand <- as.character(GenomicRanges::strand(bounds))
    result$start[pair_index] <- GenomicRanges::start(bounds)
    result$end[pair_index] <- GenomicRanges::end(bounds)
    result$orf_3p[pair_index] <- ifelse(
      path_strand == "+",
      GenomicRanges::end(bounds),
      GenomicRanges::start(bounds)
    )
    result
  }
  bounds_stop_inclusive <- pair_range_bounds(
    ranges_stop_inclusive, projected_idx, length(orf_in_tx)
  )
  bounds_stop_exclusive <- pair_range_bounds(
    ranges_stop_exclusive, projected_idx, length(orf_in_tx)
  )

  # From this point `seq_nt`, translation, and the primary genomic coordinates
  # describe the stop-inclusive chain returned in `ranges`.
  if (length(projected_idx)) {
    inclusive_start_1b <- orf_5p_tx_all[projected_idx] + 1L
    inclusive_end_1b <- inclusive_3p_tx_0b[projected_idx] + 1L
    seq_nt_chr[projected_idx] <- substring(
      tx_seq_per_orf[projected_idx],
      inclusive_start_1b,
      inclusive_end_1b
    )
    seq_aa_chr[projected_idx] <- as.character(
      Biostrings::translate(Biostrings::DNAStringSet(seq_nt_chr[projected_idx]))
    )
    genomic_st[projected_idx] <- bounds_stop_inclusive$start[projected_idx]
    genomic_en[projected_idx] <- bounds_stop_inclusive$end[projected_idx]
  }

  # --- Pre-compute lengths ---
  len_nt <- nchar(seq_nt_chr)
  len_aa <- nchar(gsub("*", "", seq_aa_chr, fixed = TRUE))

  # --- Flag out-of-frame ORF isoforms ---
  complete_codons <- len_nt %% 3 == 0
  if (any(!complete_codons)) {
    cli::cli_warn(c(
      "{sum(!complete_codons)} ORF-isoform pair{?s} {?has/have} a sequence \\
       length not divisible by 3.",
      "i" = "Flagged as {.code complete_codons = FALSE} in the output."
    ))
  }

  # --- Build status table for mapped pairs ---
  orf_status <- tibble(
    ORF_isoform_id        = ORF_isoform_id,
    ORF_id                = names(orf_by_tx),
    transcript_id         = names(tx_by_orf),
    seq_nt                = seq_nt_chr,
    seq_aa                = seq_aa_chr,
    len_nt                = len_nt,
    len_aa                = len_aa,
    seqnames              = genomic_sn,
    start                 = genomic_st,
    end                   = genomic_en,
    strand                = strands_vec,
    complete_codons       = complete_codons,
    input_stop_convention = resolved_stop_convention,
    stop_adjustment_nt    = stop_adjustment_nt,
    stop_normalization_valid = stop_normalization_valid,
    stop_normalization_status = stop_normalization_status,
    start_stop_exclusive  = bounds_stop_exclusive$start,
    end_stop_exclusive    = bounds_stop_exclusive$end,
    orf_3p_stop_exclusive = bounds_stop_exclusive$orf_3p,
    start_stop_inclusive  = bounds_stop_inclusive$start,
    end_stop_inclusive    = bounds_stop_inclusive$end,
    orf_3p_stop_inclusive = bounds_stop_inclusive$orf_3p,
    canonical_stop_codon  = canonical_stop_codon,
    canonical_complete_codons = canonical_complete_codons,
    canonical_5p_tx       = as.integer(orf_5p_tx_all + 1L),
    canonical_3p_tx       = as.integer(inclusive_3p_tx_0b + 1L),
    downstream_codon      = downstream_codon_vec,
    downstream_is_stop    = downstream_is_stop_vec,
    stop_codon_end        = stop_codon_end_vec
  ) %>%
    mutate(
      start_codon = substr(seq_nt, 1, 3),
      stop_codon  = substr(seq_nt, len_nt - 2, len_nt),
      # Status uses the canonical stop-inclusive sequence. Invalid or
      # unresolved normalisations are conservatively reported as no-stop.
      orf_status  = case_when(
        !stop_normalization_valid & !(start_codon %in% start_codons) ~
          "no_stop_no_start",
        !stop_normalization_valid ~ "no_stop",
        # Internal stop codon anywhere before the final aa
        grepl("[*]", substr(seq_aa, 1, nchar(seq_aa) - 1)) ~ "internal_stop",
        # Canonical sequence has a stop but no permitted start.
        !(start_codon %in% start_codons) ~ "no_start",
        TRUE ~ "translatable"
      )
    ) %>%
    group_by(ORF_id, orf_status) %>%
    mutate(
      unique_iso_n  = as.numeric(factor(seq_aa)),
      unique_tx_iso = paste(orf_status, unique_iso_n, sep = "_")
    ) %>%
    ungroup()

  # --- Stub rows for intronic ORF-transcript pairs ---
  if (length(intronic_idx) > 0) {
    ic               <- grl_bounds(intronic_orfs_grl)
    intronic_strands <- suppressWarnings(
      vapply(as.list(intronic_orfs_grl),
             function(x) as.character(strand(x)[1L]),
             character(1L))
    )
    intronic_stubs <- tibble(
      ORF_isoform_id     = paste(intronic_orf_ids, paste0("iso-", intronic_tx_ids), sep = "_"),
      ORF_id             = intronic_orf_ids,
      transcript_id      = intronic_tx_ids,
      seq_nt             = NA_character_,
      seq_aa             = NA_character_,
      len_nt             = NA_integer_,
      len_aa             = NA_integer_,
      seqnames           = ic$seqnames,
      start              = ic$start,
      end                = ic$end,
      strand             = intronic_strands,
      complete_codons    = NA,
      start_codon        = NA_character_,
      stop_codon         = NA_character_,
      orf_status         = "intronic",
      unique_iso_n       = NA_real_,
      unique_tx_iso      = NA_character_,
      downstream_codon   = NA_character_,
      downstream_is_stop = NA,
      stop_codon_end     = NA_integer_
    )
  } else {
    intronic_stubs <- NULL
  }

  # --- Stub rows for ORFs with no transcript overlap ---
  if (length(unmapped_orf_idx) > 0) {
    unmapped_grl  <- orfs[unmapped_orf_idx]
    uc            <- grl_bounds(unmapped_grl)
    no_tx_strands <- suppressWarnings(
      vapply(as.list(unmapped_grl),
             function(x) as.character(strand(x)[1L]),
             character(1L))
    )
    no_tx_stubs <- tibble(
      ORF_isoform_id     = paste0(names(orfs)[unmapped_orf_idx], "_iso-no_transcript"),
      ORF_id             = names(orfs)[unmapped_orf_idx],
      transcript_id      = NA_character_,
      seq_nt             = NA_character_,
      seq_aa             = NA_character_,
      len_nt             = NA_integer_,
      len_aa             = NA_integer_,
      seqnames           = uc$seqnames,
      start              = uc$start,
      end                = uc$end,
      strand             = no_tx_strands,
      complete_codons    = NA,
      start_codon        = NA_character_,
      stop_codon         = NA_character_,
      orf_status         = "no_transcript",
      unique_iso_n       = NA_real_,
      unique_tx_iso      = NA_character_,
      downstream_codon   = NA_character_,
      downstream_is_stop = NA,
      stop_codon_end     = NA_integer_
    )
  } else {
    no_tx_stubs <- NULL
  }

  # Combine all rows (mapped pairs + intronic stubs + no-transcript stubs)
  orf_status <- dplyr::bind_rows(orf_status, intronic_stubs, no_tx_stubs)

  # Global unique amino acid sequence ID.
  # Rows with the same seq_aa (across ALL ORFs and transcripts) share the same
  # unique_aa_id; stub rows with NA seq_aa are assigned NA.
  unique_seqs <- unique(orf_status$seq_aa[!is.na(orf_status$seq_aa)])
  orf_status <- orf_status %>%
    dplyr::mutate(
      unique_aa_id = dplyr::if_else(
        is.na(seq_aa),
        NA_character_,
        paste0("aa_iso_", match(seq_aa, unique_seqs))
      )
    )

  # --- Join metadata ---
  orf_status <- orf_status %>%
    left_join(transcript_meta, by = "transcript_id")

  if (!is.null(orfs_meta)) {
    orf_status <- orf_status %>%
      left_join(orfs_meta, by = "ORF_id")
  }

  # --- Optional: CDS overlap and reading-frame check ---
  # cds_gr must be a GRangesList keyed by transcript_id (same key space as
  # names(annotations$transcripts)).  Frame is computed in transcript
  # coordinates so that splice junctions between the ORF and the CDS are
  # accounted for correctly.
  if (!is.null(cds_gr)) {
    # Index from ORF_isoform_id -> position in orf_in_tx / tx_by_orf
    pair_idx <- setNames(seq_along(orf_in_tx), names(orf_in_tx))

    orf_status$overlaps_cds   <- NA
    orf_status$pct_orf_in_cds <- NA_real_
    orf_status$pct_cds_in_orf <- NA_real_
    orf_status$cds_frame      <- NA_character_
    orf_status$coding_transcript <- NA
    orf_status$cds_5p_tx <- NA_integer_
    orf_status$cds_3p_tx <- NA_integer_
    orf_status$start_in_cds_frame <- NA
    orf_status$end_in_cds_frame <- NA
    orf_status$reference_orf_type <- NA_character_
    orf_status$reference_orf_class <- NA_character_
    orf_status$reference_annotation_status <- NA_character_

    # Map each row to its pair index (NA for intronic / no-transcript stubs)
    row_pair_i <- pair_idx[orf_status$ORF_isoform_id]
    valid_rows <- which(!is.na(row_pair_i))
    orf_status$overlaps_cds[valid_rows] <- FALSE   # default for mapped rows

    tx_ids_valid <- names(tx_by_orf)[row_pair_i[valid_rows]]
    has_cds_ann  <- tx_ids_valid %in% names(cds_gr)
    if (any(has_cds_ann)) {
      has_cds_ann[has_cds_ann] <-
        S4Vectors::elementNROWS(cds_gr[tx_ids_valid[has_cds_ann]]) > 0L
    }
    orf_status$coding_transcript[valid_rows] <- has_cds_ann
    noncoding_rows <- valid_rows[!has_cds_ann]
    noncoding_pairs <- row_pair_i[noncoding_rows]
    noncoding_ready <- stop_normalization_valid[noncoding_pairs] %in% TRUE
    orf_status$reference_orf_type[noncoding_rows[noncoding_ready]] <-
      "varRNA-ORF"
    orf_status$reference_orf_class[noncoding_rows[noncoding_ready]] <-
      unname(.BILBORF_ORF_CLASS_MAP["varRNA-ORF"])
    orf_status$reference_annotation_status[noncoding_rows] <- ifelse(
      noncoding_ready,
      "compatible transcript has no annotated CDS",
      stop_normalization_status[noncoding_pairs]
    )
    valid_rows_c <- valid_rows[has_cds_ann]
    tx_ids_c     <- tx_ids_valid[has_cds_ann]
    pairs_i_c    <- row_pair_i[valid_rows_c]

    if (length(valid_rows_c) > 0L) {
      unique_tx <- unique(tx_ids_c)

      # Progress bar counts unique transcripts, not rows — much faster to tick
      pb <- utils::txtProgressBar(min = 0, max = length(unique_tx), style = 3)
      on.exit(close(pb), add = TRUE)

      for (k in seq_along(unique_tx)) {
        utils::setTxtProgressBar(pb, k)
        tid <- unique_tx[k]

        rows_k  <- valid_rows_c[tx_ids_c == tid]
        pairs_k <- pairs_i_c[tx_ids_c == tid]

        cds_exons  <- cds_gr[[tid]]
        tx_exons_k <- tx_by_orf[[pairs_k[1L]]]   # same transcript for all in group
        strnd_k    <- strands_vec[pairs_k[1L]]    # same strand for all in group
        cds_bp_k   <- sum(width(cds_exons))

        # --- Overlap & percentages: one findOverlaps for all ORFs on this tx ---
        orf_exons_list <- orf_in_tx[pairs_k]
        orf_exons_flat <- unlist(orf_exons_list, use.names = FALSE)
        mcols(orf_exons_flat) <- NULL
        orf_grp <- rep(seq_along(pairs_k), elementNROWS(orf_exons_list))

        # Exonic bp per ORF (denominator for pct_orf_in_cds)
        orf_bp_per <- as.integer(tapply(width(orf_exons_flat), orf_grp, sum))

        # Overlap with CDS exons
        ov <- findOverlaps(orf_exons_flat, cds_exons, ignore.strand = FALSE)

        olap_vec <- integer(length(pairs_k))   # 0 = no overlap
        if (length(ov) > 0L) {
          isect_w     <- width(pintersect(orf_exons_flat[queryHits(ov)],
                                          cds_exons[subjectHits(ov)]))
          olap_tbl    <- tapply(isect_w, orf_grp[queryHits(ov)], sum)
          olap_vec[as.integer(names(olap_tbl))] <- as.integer(olap_tbl)
        }

        overlapping <- olap_vec > 0L
        orf_status$overlaps_cds[rows_k[overlapping]]   <- TRUE
        orf_status$pct_orf_in_cds[rows_k[overlapping]] <-
          round(100 * olap_vec[overlapping] / orf_bp_per[overlapping], 2)
        orf_status$pct_cds_in_orf[rows_k[overlapping]] <-
          round(100 * olap_vec[overlapping] / cds_bp_k, 2)

        # --- Frame check: cds_tx_coord computed once per transcript ---
        cds_tx_coord_k <- .genomic_5prime_to_tx_coord(cds_exons, tx_exons_k, strnd_k)

        if (!is.na(cds_tx_coord_k)) {
          cds_3p_tx_coord_k <-
            cds_tx_coord_k + sum(BiocGenerics::width(cds_exons)) - 1L
          orf_5p_tx_k <- orf_5p_tx_all[pairs_k]
          orf_3p_tx_k <- inclusive_3p_tx_0b[pairs_k]

          orf_status$cds_5p_tx[rows_k] <- cds_tx_coord_k + 1L
          orf_status$cds_3p_tx[rows_k] <- cds_3p_tx_coord_k + 1L

          # Batch ORF 5' genomic positions using tapply on the already-flat vector
          g_vec_k <- if (strnd_k == "+") {
            as.integer(tapply(start(orf_exons_flat), orf_grp, min))
          } else {
            as.integer(tapply(end(orf_exons_flat),   orf_grp, max))
          }

          # Vectorised transcript-coordinate lookup (findInterval-based)
          orf_tx_coords_k <- .genomic_5prime_to_tx_coord_vec(g_vec_k, tx_exons_k, strnd_k)

          valid_frame <- !is.na(orf_tx_coords_k)
          orf_status$cds_frame[rows_k[valid_frame]] <- ifelse(
            (orf_tx_coords_k[valid_frame] - cds_tx_coord_k) %% 3L == 0L,
            "in_frame", "out_of_frame"
          )

          classification_ready <-
            stop_normalization_valid[pairs_k] %in% TRUE &
            !is.na(orf_5p_tx_k) & !is.na(orf_3p_tx_k)
          if (any(classification_ready)) {
            classified_rows <- rows_k[classification_ready]
            classified_types <- .classify_orf_cds_geometry(
              orf_5p_tx_k[classification_ready],
              orf_3p_tx_k[classification_ready],
              rep(cds_tx_coord_k, sum(classification_ready)),
              rep(cds_3p_tx_coord_k, sum(classification_ready))
            )
            start_frame <-
              (orf_5p_tx_k[classification_ready] - cds_tx_coord_k) %% 3L == 0L
            end_frame <-
              (orf_3p_tx_k[classification_ready] - cds_3p_tx_coord_k) %% 3L == 0L

            orf_status$start_in_cds_frame[classified_rows] <- start_frame
            orf_status$end_in_cds_frame[classified_rows] <- end_frame
            orf_status$reference_orf_type[classified_rows] <- classified_types
            orf_status$reference_orf_class[classified_rows] <- unname(
              .BILBORF_ORF_CLASS_MAP[classified_types]
            )
            orf_status$reference_annotation_status[classified_rows] <-
              "classified against annotated CDS"
          }

          not_ready <- !classification_ready
          orf_status$reference_annotation_status[rows_k[not_ready]] <-
            stop_normalization_status[pairs_k[not_ready]]
        } else {
          orf_status$reference_annotation_status[rows_k] <-
            "CDS could not be projected onto transcript"
        }
      }

      close(pb)
      on.exit()   # clear safety handler after explicit close
    }
  }

  # --- Column order ---
  orf_status <- orf_status %>%
    dplyr::relocate(ORF_isoform_id, ORF_id, transcript_id, gene_id, gene_name,
                    seqnames, start, end, strand, start_codon, stop_codon,
                    orf_status, unique_aa_id, complete_codons,
                    dplyr::any_of(c(
                      "input_stop_convention", "stop_adjustment_nt",
                      "stop_normalization_valid", "stop_normalization_status",
                      "start_stop_exclusive", "end_stop_exclusive",
                      "orf_3p_stop_exclusive", "start_stop_inclusive",
                      "end_stop_inclusive", "orf_3p_stop_inclusive",
                      "canonical_stop_codon", "canonical_complete_codons",
                      "canonical_5p_tx", "canonical_3p_tx"
                    )),
                    dplyr::any_of(c("downstream_codon", "downstream_is_stop",
                                    "stop_codon_end", "overlaps_cds",
                                    "pct_orf_in_cds", "pct_cds_in_orf",
                                    "cds_frame", "coding_transcript",
                                    "cds_5p_tx", "cds_3p_tx",
                                    "start_in_cds_frame", "end_in_cds_frame",
                                    "reference_orf_type",
                                    "reference_orf_class",
                                    "reference_annotation_status")))

  if ("reference_orf_type" %in% names(orf_status)) {
    orf_status$reference_orf_type <- factor(
      orf_status$reference_orf_type,
      levels = .BILBORF_ORF_TYPE_LEVELS
    )
    orf_status$reference_orf_class <- factor(
      orf_status$reference_orf_class,
      levels = c("canonical", "variant of canonical", "non-canonical")
    )
  }

  return(list(table = orf_status, ranges = ranges_stop_inclusive))
}
