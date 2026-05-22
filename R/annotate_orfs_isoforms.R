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
#' @param check_stop_codon Logical. If \code{TRUE}, checks whether the 3 nt
#'   immediately 3' of each ORF on its mapped transcript form a stop codon.
#'   Splice junctions are handled correctly: if the ORF ends in one exon the
#'   search continues into the next exon. Adds three columns:
#'   \code{downstream_codon} (the 3 nt, or \code{NA} if the transcript ends
#'   within 3 nt), \code{downstream_is_stop} (logical), and
#'   \code{stop_codon_end} (strand-aware genomic position of the last nt of
#'   the stop codon: \code{max(end)} for + strand, \code{min(start)} for -
#'   strand). Default \code{FALSE}.
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
#'   \code{NA}). Default \code{NULL}.
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
#'     \item{ranges}{A GRangesList with genomic coordinates of ORF-transcript
#'       exonic intersections. Contains entries only for rows with actual
#'       sequence (\code{orf_status} is not \code{"intronic"} or
#'       \code{"no_transcript"}).}
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
#' @note ORF coordinates should represent the complete ORF boundaries (start to
#'   stop codon), not just the coding sequence. For multi-exon ORFs, provide
#'   all exonic segments.
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
                                  check_stop_codon = FALSE,
                                  cds_gr           = NULL,
                                  chunk_size       = NULL) {

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

  # Helper: min/max genomic bounds per element of a GRangesList
  grl_bounds <- function(grl) {
    ul  <- unlist(grl, use.names = FALSE)
    grp <- rep(seq_along(grl), elementNROWS(grl))
    fst <- cumsum(c(1L, elementNROWS(grl)[-length(grl)]))
    list(
      seqnames = as.character(seqnames(ul)[fst]),
      start    = as.integer(tapply(start(ul), grp, min)),
      end      = as.integer(tapply(end(ul),   grp, max))
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
  seq_aa     <- suppressWarnings(translate(seq_nt))
  seq_nt_chr <- as.character(seq_nt)
  seq_aa_chr <- as.character(seq_aa)

  # --- Vectorised genomic coordinates for mapped pairs ---
  bc         <- grl_bounds(orf_in_tx)
  genomic_sn <- bc$seqnames
  genomic_st <- bc$start
  genomic_en <- bc$end

  # Strand for each ORF-tx pair (all exons share the same strand)
  strands_vec <- suppressWarnings(
    vapply(seq_len(length(orf_in_tx)),
           function(i) as.character(strand(orf_in_tx[[i]])[1L]),
           character(1L))
  )

  # --- Optional: downstream stop codon check ---
  downstream_codon_vec   <- rep(NA_character_, length(orf_in_tx))
  downstream_is_stop_vec <- rep(NA,             length(orf_in_tx))
  stop_codon_end_vec     <- rep(NA_integer_,    length(orf_in_tx))

  if (check_stop_codon) {
    ds_gr_list <- suppressWarnings(
      mapply(.get_downstream_nt,
             as.list(orf_in_tx), as.list(tx_by_orf),
             MoreArgs = list(n = 3L), SIMPLIFY = FALSE)
    )
    has_ds <- lengths(ds_gr_list) > 0L
    if (any(has_ds)) {
      ds_grl  <- GRangesList(ds_gr_list[has_ds])
      ds_seqs <- as.character(extractTranscriptSeqs(BSgenome, ds_grl))
      downstream_codon_vec[has_ds]   <- ds_seqs
      downstream_is_stop_vec[has_ds] <- ds_seqs %in% stop_codons
      for (i in which(has_ds)) {
        gr <- ds_gr_list[[i]]
        stop_codon_end_vec[i] <-
          if (suppressWarnings(as.character(strand(gr)[1L])) == "+") max(end(gr)) else min(start(gr))
      }
    }
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
    downstream_codon      = downstream_codon_vec,
    downstream_is_stop    = downstream_is_stop_vec,
    stop_codon_end        = stop_codon_end_vec
  ) %>%
    mutate(
      start_codon = substr(seq_nt, 1, 3),
      stop_codon  = substr(seq_nt, len_nt - 2, len_nt),
      # orf_status is assigned after downstream_is_stop is available so that
      # an ORF whose stop codon lies immediately downstream (not included in its
      # coordinates) is still classified correctly.
      # !(downstream_is_stop %in% TRUE) is TRUE for both FALSE and NA, so when
      # check_stop_codon = FALSE (all NA) the logic is identical to before.
      orf_status  = case_when(
        # Internal stop codon anywhere before the final aa
        grepl("[*]", substr(seq_aa, 1, nchar(seq_aa) - 1)) ~ "internal_stop",
        # No start AND no stop (in seq or downstream)
        !(start_codon %in% start_codons) & !(stop_codon %in% stop_codons) &
          !(downstream_is_stop %in% TRUE) ~ "no_stop_no_start",
        # Has start but no stop (in seq or downstream)
        !(stop_codon  %in% stop_codons) &
          !(downstream_is_stop %in% TRUE) ~ "no_stop",
        # Has stop (in seq or downstream) but no start
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

    orf_status$overlaps_cds    <- NA
    orf_status$pct_orf_in_cds  <- NA_real_
    orf_status$pct_cds_in_orf  <- NA_real_
    orf_status$cds_frame       <- NA_character_

    n_rows  <- nrow(orf_status)
    pb      <- utils::txtProgressBar(min = 0, max = n_rows, style = 3)
    on.exit(close(pb), add = TRUE)
    for (row in seq_len(n_rows)) {
      utils::setTxtProgressBar(pb, row)
      iid <- orf_status$ORF_isoform_id[row]
      i   <- pair_idx[iid]
      if (is.na(i)) next            # intronic or no-transcript stub

      tx_id <- names(tx_by_orf)[i]
      orf_status$overlaps_cds[row] <- FALSE   # default once transcript found
      if (!tx_id %in% names(cds_gr)) next

      cds_exons   <- cds_gr[[tx_id]]
      orf_exons_i <- orf_in_tx[[i]]
      tx_exons_i  <- tx_by_orf[[i]]
      strnd       <- strands_vec[i]

      # Exon-level overlap; pintersect gives bp overlap per hit
      ov <- findOverlaps(orf_exons_i, cds_exons, ignore.strand = FALSE)
      if (length(ov) == 0L) next

      orf_status$overlaps_cds[row] <- TRUE

      # Percentage overlap (splice-aware: use exonic bp, not genomic span)
      overlap_bp <- sum(width(pintersect(
        orf_exons_i[queryHits(ov)],
        cds_exons[subjectHits(ov)]
      )))
      orf_bp <- sum(width(orf_exons_i))
      cds_bp <- sum(width(cds_exons))
      orf_status$pct_orf_in_cds[row] <- round(100 * overlap_bp / orf_bp, 2)
      orf_status$pct_cds_in_orf[row] <- round(100 * overlap_bp / cds_bp, 2)

      # Transcript-coordinate positions of the two 5' ends
      orf_tx_coord <- .genomic_5prime_to_tx_coord(orf_exons_i, tx_exons_i, strnd)
      cds_tx_coord <- .genomic_5prime_to_tx_coord(cds_exons,   tx_exons_i, strnd)

      if (!is.na(orf_tx_coord) && !is.na(cds_tx_coord)) {
        orf_status$cds_frame[row] <-
          if ((orf_tx_coord - cds_tx_coord) %% 3L == 0L) "in_frame" else "out_of_frame"
      }
    }
    close(pb)
    on.exit()    # clear the safety on.exit now we closed it explicitly
  }

  # --- Column order ---
  orf_status <- orf_status %>%
    dplyr::relocate(ORF_isoform_id, ORF_id, transcript_id, gene_id, gene_name,
                    seqnames, start, end, strand, start_codon, stop_codon,
                    orf_status, unique_aa_id, complete_codons,
                    dplyr::any_of(c("downstream_codon", "downstream_is_stop",
                                    "stop_codon_end", "overlaps_cds",
                                    "pct_orf_in_cds", "pct_cds_in_orf",
                                    "cds_frame")))

  return(list(table = orf_status, ranges = orf_in_tx))
}
