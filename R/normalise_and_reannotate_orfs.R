#' Normalise ORF stops and re-annotate isoforms in one pass
#'
#' Reuse the ORF-transcript mappings produced by [reannotate_orf_type()] to
#' validate an included or omitted stop codon, select a reference isoform, and
#' return the corresponding canonical stop-inclusive exon chain. This avoids
#' rebuilding exon intersections separately for stop normalisation and type
#' annotation.
#'
#' @param orfs A `GRanges` of unique ORF spans or a `GRangesList` of unique ORF
#'   chains. Coordinates are 1-based and inclusive.
#' @param transcripts A named `GRangesList` of transcript exons.
#' @param cds_by_tx A named `GRangesList` of CDS exons.
#' @param genome A genome accepted by
#'   [GenomicFeatures::extractTranscriptSeqs()], such as a `BSgenome`.
#' @param id_col Optional ORF-ID metadata column; otherwise `names(orfs)` are
#'   used.
#' @param type_priority Priority used to select among valid transcript
#'   interpretations.
#' @param stop_codon_convention Whether the supplied 3' boundary includes or
#'   excludes the stop codon. `"auto"` accepts a convention only when exactly
#'   one of the endpoint triplet and following triplet is a canonical stop.
#' @param stop_codons Accepted stop triplets in transcript orientation.
#' @param allow_terminal_extrapolation,max_terminal_overhang Passed to
#'   [reannotate_orf_type()]. Stop-codon normalisation itself requires the
#'   canonical endpoint to exist on a compatible transcript; a longer isoform
#'   is therefore preferred over extrapolating an absent terminal exon.
#' @param unmatched_type Passed to [reannotate_orf_type()].
#'
#' @return A list with `orfs`, `pairs`, and `canonical_chains`. `orfs` contains
#'   one row per input ORF. `pairs` retains every candidate ORF-transcript
#'   interpretation and its stop-normalisation audit. `canonical_chains` is a
#'   stop-inclusive `GRangesList` for successfully normalised ORFs, named by
#'   ORF ID.
#' @export
normalise_and_reannotate_orfs <- function(
    orfs,
    transcripts,
    cds_by_tx,
    genome,
    id_col = NULL,
    type_priority = .BILBORF_REFERENCE_ORF_TYPE_PRIORITY,
    stop_codon_convention = c("excluded", "included", "auto"),
    stop_codons = c("TAA", "TAG", "TGA"),
    allow_terminal_extrapolation = TRUE,
    max_terminal_overhang = 30L,
    unmatched_type = NA_character_) {

  stop_codon_convention <- match.arg(stop_codon_convention)
  stop_codons <- toupper(as.character(stop_codons))
  if (!length(stop_codons) || anyNA(stop_codons) ||
      any(nchar(stop_codons) != 3L)) {
    stop("stop_codons must contain non-missing three-base strings",
         call. = FALSE)
  }

  annotation <- reannotate_orf_type(
    orfs = orfs,
    transcripts = transcripts,
    cds_by_tx = cds_by_tx,
    id_col = id_col,
    type_priority = type_priority,
    stop_codon_convention = stop_codon_convention,
    allow_terminal_extrapolation = allow_terminal_extrapolation,
    max_terminal_overhang = max_terminal_overhang,
    unmatched_type = unmatched_type
  )

  pairs <- annotation$pairs
  if (!nrow(pairs)) {
    annotation$pairs$canonical_5p_tx <- integer()
    annotation$pairs$canonical_3p_tx <- integer()
    annotation$pairs$canonical_stop_triplet <- character()
    annotation$pairs$canonical_stop_valid <- logical()
    annotation$pairs$canonical_complete_codons <- logical()
    annotation$pairs$stop_normalization_status <- character()
    annotation$canonical_chains <- GenomicRanges::GRangesList()
    return(annotation)
  }

  used_tx <- unique(pairs$transcript_id)
  tx_sequences <- GenomicFeatures::extractTranscriptSeqs(
    genome, transcripts[used_tx]
  )
  names(tx_sequences) <- used_tx
  tx_width <- stats::setNames(nchar(as.character(tx_sequences)), used_tx)

  raw_3p <- as.integer(pairs$orf_3p_tx_effective)
  raw_5p <- as.integer(pairs$orf_5p_tx_effective)
  current_triplet <- .pair_transcript_triplet(
    pairs$transcript_id, raw_3p, tx_sequences, offset = 0L
  )
  following_triplet <- .pair_transcript_triplet(
    pairs$transcript_id, raw_3p, tx_sequences, offset = 3L
  )
  current_valid <- current_triplet %in% stop_codons
  following_valid <- following_triplet %in% stop_codons

  adjustment <- rep(NA_integer_, nrow(pairs))
  normalization_status <- rep("stop convention unresolved", nrow(pairs))
  if (stop_codon_convention == "included") {
    adjustment[] <- 0L
    normalization_status[] <- "input stop included"
  } else if (stop_codon_convention == "excluded") {
    adjustment[] <- 3L
    normalization_status[] <- "input stop excluded; added 3 transcript nt"
  } else {
    included_only <- current_valid & !following_valid
    excluded_only <- following_valid & !current_valid
    adjustment[included_only] <- 0L
    adjustment[excluded_only] <- 3L
    normalization_status[included_only] <- "auto: input stop included"
    normalization_status[excluded_only] <-
      "auto: input stop excluded; added 3 transcript nt"
    normalization_status[current_valid & following_valid] <-
      "auto: both stop conventions are plausible"
    normalization_status[!current_valid & !following_valid] <-
      "auto: neither stop convention has a canonical stop"
  }

  canonical_3p <- raw_3p + adjustment
  canonical_stop_triplet <- ifelse(
    adjustment == 0L, current_triplet,
    ifelse(adjustment == 3L, following_triplet, NA_character_)
  )
  interval_valid <-
    !is.na(raw_5p) & !is.na(canonical_3p) &
    raw_5p >= 1L & canonical_3p >= raw_5p &
    canonical_3p <= unname(tx_width[pairs$transcript_id])
  canonical_width <- canonical_3p - raw_5p + 1L
  complete_codons <- interval_valid & canonical_width %% 3L == 0L
  valid_stop <- interval_valid & canonical_stop_triplet %in% stop_codons
  normalization_valid <- valid_stop & complete_codons

  normalization_status[!is.na(adjustment) & !interval_valid] <-
    "canonical stop lies outside the compatible transcript"
  normalization_status[interval_valid & !valid_stop] <-
    "canonical endpoint is not a permitted stop codon"
  normalization_status[valid_stop & !complete_codons] <-
    "canonical ORF width is not divisible by three"

  pairs$canonical_5p_tx <- raw_5p
  pairs$canonical_3p_tx <- as.integer(canonical_3p)
  pairs$canonical_stop_triplet <- canonical_stop_triplet
  pairs$canonical_stop_valid <- valid_stop
  pairs$canonical_width_nt <- as.integer(canonical_width)
  pairs$canonical_complete_codons <- complete_codons
  pairs$stop_normalization_valid <- normalization_valid
  pairs$stop_normalization_status <- normalization_status

  # Recompute coding classifications from the validated canonical endpoint.
  # This is essential for sequence-driven `auto` normalisation and also makes
  # the order explicit: normalise first, classify second, using the same pair
  # mapping produced above.
  classification_ready <-
    normalization_valid &
    pairs$splice_chain_compatible &
    (pairs$boundaries_exonic | pairs$terminal_extrapolation_used) &
    pairs$coding_transcript &
    !is.na(pairs$cds_5p_tx) & !is.na(pairs$cds_3p_tx)
  if (any(classification_ready)) {
    pairs$reference_orf_type[classification_ready] <-
      .classify_orf_cds_geometry(
        pairs$canonical_5p_tx[classification_ready],
        pairs$canonical_3p_tx[classification_ready],
        pairs$cds_5p_tx[classification_ready],
        pairs$cds_3p_tx[classification_ready]
      )
    pairs$start_in_cds_frame[classification_ready] <-
      (pairs$canonical_5p_tx[classification_ready] -
       pairs$cds_5p_tx[classification_ready]) %% 3L == 0L
    pairs$end_in_cds_frame[classification_ready] <-
      (pairs$canonical_3p_tx[classification_ready] -
       pairs$cds_3p_tx[classification_ready]) %% 3L == 0L
  }

  valid_pairs <- pairs |>
    dplyr::filter(
      splice_chain_compatible,
      stop_normalization_valid,
      !is.na(reference_orf_type)
    ) |>
    dplyr::group_by(orf_id) |>
    dplyr::filter(
      !any(!terminal_extrapolation_used) | !terminal_extrapolation_used
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(.priority = match(reference_orf_type, type_priority))

  best_pair <- valid_pairs |>
    dplyr::arrange(orf_id, .priority, transcript_id) |>
    dplyr::group_by(orf_id) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup()

  selected <- best_pair |>
    dplyr::transmute(
      orf_id,
      matched_transcript_id = transcript_id,
      reference_orf_type,
      reference_start_in_frame = start_in_cds_frame,
      reference_end_in_frame = end_in_cds_frame,
      reference_splice_chain_checked = splice_chain_checked,
      reference_splice_chain_compatible = splice_chain_compatible,
      reference_boundary_extrapolated = terminal_extrapolation_used,
      reference_5p_overhang_nt = overhang_5p_nt,
      reference_3p_overhang_nt = overhang_3p_nt,
      reference_stop_boundary_adjusted = canonical_3p_tx != orf_3p_tx_effective,
      reference_stop_boundary_adjustment_nt =
        canonical_3p_tx - orf_3p_tx_effective,
      reference_annotation_status = annotation_status,
      canonical_5p_tx,
      canonical_3p_tx,
      canonical_stop_triplet,
      canonical_stop_valid,
      canonical_width_nt,
      canonical_complete_codons,
      stop_normalization_status
    )

  pair_summary <- valid_pairs |>
    dplyr::arrange(orf_id, .priority) |>
    dplyr::group_by(orf_id) |>
    dplyr::summarise(
      n_compatible_transcripts = dplyr::n_distinct(transcript_id),
      n_reference_orf_types = dplyr::n_distinct(reference_orf_type),
      all_reference_orf_types = paste(unique(reference_orf_type), collapse = ","),
      reference_type_ambiguous = n_reference_orf_types > 1L,
      .groups = "drop"
    )

  failure_status <- pairs |>
    dplyr::group_by(orf_id) |>
    dplyr::summarise(
      stop_failure = dplyr::case_when(
        any(stop_normalization_valid & splice_chain_compatible) ~ NA_character_,
        any(splice_chain_compatible &
              stop_normalization_status ==
                "canonical stop lies outside the compatible transcript") ~
          "no compatible transcript contains the complete stop codon",
        any(splice_chain_compatible) ~
          "no compatible transcript has a valid canonical stop",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    )

  base_orfs <- annotation$orfs |>
    dplyr::select(
      orf_id, seqnames, start, end, strand,
      n_orf_exons, splice_chain_supplied
    )
  original_status <- stats::setNames(
    annotation$orfs$reference_annotation_status,
    annotation$orfs$orf_id
  )
  orf_table <- base_orfs |>
    dplyr::left_join(selected, by = "orf_id") |>
    dplyr::left_join(pair_summary, by = "orf_id") |>
    dplyr::left_join(failure_status, by = "orf_id") |>
    dplyr::mutate(
      reference_annotation_status = dplyr::coalesce(
        reference_annotation_status,
        stop_failure,
        unname(original_status[orf_id])
      ),
      reference_boundary_extrapolated = dplyr::coalesce(
        reference_boundary_extrapolated, FALSE
      ),
      reference_stop_boundary_adjusted = dplyr::coalesce(
        reference_stop_boundary_adjusted, FALSE
      ),
      reference_stop_boundary_adjustment_nt = dplyr::coalesce(
        reference_stop_boundary_adjustment_nt, 0L
      ),
      canonical_stop_valid = dplyr::coalesce(canonical_stop_valid, FALSE),
      canonical_complete_codons = dplyr::coalesce(
        canonical_complete_codons, FALSE
      ),
      n_compatible_transcripts = dplyr::coalesce(
        n_compatible_transcripts, 0L
      ),
      n_reference_orf_types = dplyr::coalesce(n_reference_orf_types, 0L),
      reference_type_ambiguous = dplyr::coalesce(
        reference_type_ambiguous, FALSE
      ),
      reference_orf_type = factor(
        reference_orf_type, levels = .BILBORF_ORF_TYPE_LEVELS
      ),
      reference_orf_class = factor(
        unname(.BILBORF_ORF_CLASS_MAP[as.character(reference_orf_type)]),
        levels = c("canonical", "variant of canonical", "non-canonical")
      )
    ) |>
    dplyr::select(-stop_failure)

  canonical_paths <- lapply(seq_len(nrow(best_pair)), function(i) {
    .slice_transcript_interval(
      transcripts[[best_pair$transcript_id[i]]],
      best_pair$canonical_5p_tx[i],
      best_pair$canonical_3p_tx[i]
    )
  })
  canonical_chains <- GenomicRanges::GRangesList(canonical_paths)
  names(canonical_chains) <- best_pair$orf_id

  list(
    orfs = orf_table,
    pairs = pairs,
    canonical_chains = canonical_chains
  )
}

.pair_transcript_triplet <- function(
    transcript_id, endpoint, transcript_sequences, offset = 0L) {
  answer <- rep(NA_character_, length(endpoint))
  groups <- split(seq_along(endpoint), transcript_id)
  for (tx_id in names(groups)) {
    idx <- groups[[tx_id]]
    sequence <- as.character(transcript_sequences[[tx_id]])
    triplet_end <- endpoint[idx] + offset
    valid <- !is.na(triplet_end) & triplet_end >= 3L &
      triplet_end <= nchar(sequence)
    answer[idx[valid]] <- substring(
      sequence,
      triplet_end[valid] - 2L,
      triplet_end[valid]
    )
  }
  toupper(answer)
}

.slice_transcript_interval <- function(transcript, from, to) {
  from <- as.integer(from)
  to <- as.integer(to)
  if (length(from) != 1L || length(to) != 1L ||
      is.na(from) || is.na(to) || from < 1L || to < from) {
    stop("Invalid transcript interval", call. = FALSE)
  }

  tx <- ORFik::sortPerGroup(
    GenomicRanges::GRangesList(tx = transcript)
  )[[1L]]
  exon_width <- BiocGenerics::width(tx)
  exon_tx_end <- cumsum(exon_width)
  exon_tx_start <- exon_tx_end - exon_width + 1L
  keep <- exon_tx_end >= from & exon_tx_start <= to
  if (!any(keep) || to > sum(exon_width)) {
    stop("Transcript interval lies outside the transcript", call. = FALSE)
  }

  tx <- tx[keep]
  local_from <- pmax(from, exon_tx_start[keep]) - exon_tx_start[keep] + 1L
  local_to <- pmin(to, exon_tx_end[keep]) - exon_tx_start[keep] + 1L
  tx_strand <- unique(as.character(GenomicRanges::strand(tx)))
  if (length(tx_strand) != 1L || !tx_strand %in% c("+", "-")) {
    stop("Transcript must have one explicit strand", call. = FALSE)
  }

  if (tx_strand == "+") {
    genomic_start <- GenomicRanges::start(tx) + local_from - 1L
    genomic_end <- GenomicRanges::start(tx) + local_to - 1L
  } else {
    genomic_start <- GenomicRanges::end(tx) - local_to + 1L
    genomic_end <- GenomicRanges::end(tx) - local_from + 1L
  }

  path <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(tx),
    ranges = IRanges::IRanges(genomic_start, genomic_end),
    strand = GenomicRanges::strand(tx)
  )
  S4Vectors::mcols(path) <- NULL
  ORFik::sortPerGroup(GenomicRanges::GRangesList(path = path))[[1L]]
}
