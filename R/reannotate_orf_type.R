# Accepted output vocabulary. Kept internal so the package API can expose a
# copy through orf_type_levels() without allowing callers to mutate it.
.BILBORF_ORF_TYPE_LEVELS <- c(
  "annotated CDS", "uORF", "uoORF", "dORF", "doORF", "intORF",
  "N-terminal extension", "N-terminal truncation",
  "C-terminal extension", "C-terminal truncation",
  "NC-terminal extension", "varRNA-ORF"
)

.BILBORF_REFERENCE_ORF_TYPE_PRIORITY <- c(
  "annotated CDS",
  "N-terminal extension", "N-terminal truncation",
  "C-terminal extension", "C-terminal truncation",
  "NC-terminal extension",
  "intORF", "uoORF", "doORF", "uORF", "dORF",
  "varRNA-ORF"
)

.BILBORF_ORF_CLASS_MAP <- c(
  "annotated CDS" = "canonical",
  "N-terminal extension" = "variant of canonical",
  "N-terminal truncation" = "variant of canonical",
  "C-terminal extension" = "variant of canonical",
  "C-terminal truncation" = "variant of canonical",
  "NC-terminal extension" = "variant of canonical",
  "uORF" = "non-canonical",
  "uoORF" = "non-canonical",
  "dORF" = "non-canonical",
  "doORF" = "non-canonical",
  "intORF" = "non-canonical",
  "varRNA-ORF" = "non-canonical"
)

#' Supported reference ORF types
#'
#' Return the ordered vocabulary used by [reannotate_orf_type()].
#'
#' @return A character vector of supported ORF types.
#' @export
orf_type_levels <- function() {
  .BILBORF_ORF_TYPE_LEVELS
}

.classify_orf_cds_geometry <- function(orf_5p, orf_3p, cds_5p, cds_3p) {
  start_in_frame <- (orf_5p - cds_5p) %% 3L == 0L
  end_in_frame <- (orf_3p - cds_3p) %% 3L == 0L

  dplyr::case_when(
    orf_5p == cds_5p & orf_3p == cds_3p ~ "annotated CDS",
    orf_5p < cds_5p & orf_3p == cds_3p & start_in_frame ~
      "N-terminal extension",
    orf_5p > cds_5p & orf_5p <= cds_3p & orf_3p == cds_3p &
      start_in_frame ~ "N-terminal truncation",
    orf_5p == cds_5p & orf_3p > cds_3p & end_in_frame ~
      "C-terminal extension",
    orf_5p == cds_5p & orf_3p < cds_3p & orf_3p >= cds_5p &
      end_in_frame ~ "C-terminal truncation",
    orf_5p < cds_5p & orf_3p > cds_3p & start_in_frame & end_in_frame ~
      "NC-terminal extension",
    orf_3p < cds_5p ~ "uORF",
    orf_5p < cds_5p & orf_3p >= cds_5p ~ "uoORF",
    orf_5p > cds_3p ~ "dORF",
    orf_5p <= cds_3p & orf_3p > cds_3p ~ "doORF",
    orf_5p >= cds_5p & orf_3p <= cds_3p ~ "intORF",
    .default = "varRNA-ORF"
  )
}

.validate_named_grl <- function(x, arg) {
  if (!methods::is(x, "GRangesList")) {
    stop(arg, " must be a GRangesList", call. = FALSE)
  }
  if (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x))) {
    stop(arg, " must have unique, non-missing transcript-ID names",
         call. = FALSE)
  }
  invisible(x)
}

#' Re-annotate ORF type from transcript and CDS geometry
#'
#' Classify genomic ORFs relative to every compatible reference transcript and
#' its annotated CDS. Genomic boundaries are mapped to spliced transcript
#' coordinates with ORFik, so exon junctions and negative-strand transcripts
#' are handled in transcript orientation.
#'
#' @param orfs A `GRanges` with one genomic bounding range per ORF. Each range
#'   must have an explicit `+` or `-` strand. Coordinates must use the same
#'   inclusive boundary convention as `cds_by_tx`.
#' @param transcripts A `GRangesList` of transcript exons with unique transcript
#'   IDs in `names(transcripts)`.
#' @param cds_by_tx A `GRangesList` of CDS exons named in the same transcript-ID
#'   space as `transcripts`. Non-coding transcripts may be absent.
#' @param id_col Optional metadata-column name in `mcols(orfs)` containing
#'   unique ORF IDs. When `NULL`, `names(orfs)` are used; sequential IDs are
#'   generated only when names are absent.
#' @param type_priority Character vector containing every value returned by
#'   [orf_type_levels()] exactly once. It determines which interpretation wins
#'   when an ORF has different types on different transcript isoforms.
#'
#' @return A list with two tibbles:
#' \describe{
#'   \item{`orfs`}{One row per input ORF, including the selected
#'     `reference_orf_type`, `reference_orf_class`, selected transcript,
#'     frame flags, and ambiguity summaries.}
#'   \item{`pairs`}{One row per candidate ORF-transcript pair, retaining all
#'     projected coordinates, classifications, and failure statuses for audit.}
#' }
#'
#' @details
#' A transcript is compatible when both ORF boundaries map to its exons. The
#' following mutually exclusive rules are evaluated in transcript coordinates:
#'
#' * exact CDS boundaries: `annotated CDS`;
#' * same CDS 3' boundary and an in-frame upstream/downstream 5' boundary:
#'   `N-terminal extension` or `N-terminal truncation`;
#' * same CDS 5' boundary and an in-frame downstream/upstream 3' boundary:
#'   `C-terminal extension` or `C-terminal truncation`;
#' * in-frame extension beyond both CDS boundaries: `NC-terminal extension`;
#' * wholly upstream/downstream: `uORF`/`dORF`;
#' * overlapping the CDS from upstream/downstream: `uoORF`/`doORF`;
#' * contained within the CDS after the exact/truncation rules: `intORF`;
#' * a compatible transcript without a CDS, or no compatible transcript:
#'   `varRNA-ORF`.
#'
#' Exact CDS and in-frame canonical variants take precedence over a
#' non-canonical interpretation on another isoform by default. This function
#' uses geometry only; it does not inspect genomic sequence for start codons,
#' stop codons, or internal stops.
#'
#' @examples
#' library(GenomicRanges)
#' tx <- GRangesList(
#'   tx1 = GRanges("chr1", IRanges(100, 499), strand = "+")
#' )
#' cds <- GRangesList(
#'   tx1 = GRanges("chr1", IRanges(200, 399), strand = "+")
#' )
#' orfs <- GRanges("chr1", IRanges(c(200, 197), c(399, 399)), strand = "+")
#' names(orfs) <- c("exact", "n_extension")
#'
#' annotation <- reannotate_orf_type(orfs, tx, cds)
#' annotation$orfs[, c("orf_id", "reference_orf_type")]
#'
#' @export
reannotate_orf_type <- function(
    orfs,
    transcripts,
    cds_by_tx,
    id_col = NULL,
    type_priority = .BILBORF_REFERENCE_ORF_TYPE_PRIORITY) {

  if (!methods::is(orfs, "GRanges")) {
    stop("orfs must be a GRanges", call. = FALSE)
  }
  .validate_named_grl(transcripts, "transcripts")
  .validate_named_grl(cds_by_tx, "cds_by_tx")

  if (!setequal(type_priority, .BILBORF_ORF_TYPE_LEVELS) ||
      anyDuplicated(type_priority) ||
      length(type_priority) != length(.BILBORF_ORF_TYPE_LEVELS)) {
    stop("type_priority must contain every supported ORF type exactly once",
         call. = FALSE)
  }
  if (any(!as.character(GenomicRanges::strand(orfs)) %in% c("+", "-"))) {
    stop("Every ORF must have an explicit '+' or '-' strand", call. = FALSE)
  }

  orf_ids <- if (!is.null(id_col)) {
    if (!id_col %in% colnames(GenomicRanges::mcols(orfs))) {
      stop("id_col is not present in mcols(orfs): ", id_col, call. = FALSE)
    }
    as.character(GenomicRanges::mcols(orfs)[[id_col]])
  } else if (!is.null(names(orfs)) && all(nzchar(names(orfs)))) {
    names(orfs)
  } else {
    paste0("ORF_", seq_along(orfs))
  }
  if (anyNA(orf_ids) || any(!nzchar(orf_ids)) || anyDuplicated(orf_ids)) {
    stop("ORF IDs must be unique, non-missing and non-empty", call. = FALSE)
  }

  tx_ids <- names(transcripts)
  tx_spans <- unlist(range(transcripts), use.names = FALSE)
  names(tx_spans) <- tx_ids

  hits <- GenomicRanges::findOverlaps(
    orfs, tx_spans, ignore.strand = FALSE
  )
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)

  pair_template <- tibble::tibble(
    orf_id = character(), transcript_id = character(),
    coding_transcript = logical(), boundaries_exonic = logical(),
    orf_5p_tx = integer(), orf_3p_tx = integer(),
    cds_5p_tx = integer(), cds_3p_tx = integer(),
    start_in_cds_frame = logical(), end_in_cds_frame = logical(),
    reference_orf_type = character(), annotation_status = character()
  )

  if (!length(qh)) {
    pair_table <- pair_template
  } else {
    pair_orfs <- orfs[qh]
    pair_tx <- transcripts[sh]
    pair_tx_ids <- tx_ids[sh]

    # ORFik 1.29.x requires grouped ranges for startSites()/stopSites().
    # Each candidate pair is deliberately its own one-range group.
    pair_orfs_grl <- ORFik::groupGRangesBy(
      pair_orfs, seq_along(pair_orfs)
    )

    orf_5p_g <- ORFik::startSites(
      pair_orfs_grl, asGR = TRUE, keep.names = FALSE, is.sorted = TRUE
    )
    orf_3p_g <- ORFik::stopSites(
      pair_orfs_grl, asGR = TRUE, keep.names = FALSE, is.sorted = TRUE
    )
    orf_5p_mapped <- ORFik::pmapToTranscriptF(
      orf_5p_g, pair_tx, x.is.sorted = TRUE, tx.is.sorted = FALSE,
      set.seqlengths = FALSE
    )
    orf_3p_mapped <- ORFik::pmapToTranscriptF(
      orf_3p_g, pair_tx, x.is.sorted = TRUE, tx.is.sorted = FALSE,
      set.seqlengths = FALSE
    )

    orf_5p_tx <- as.integer(GenomicRanges::start(orf_5p_mapped))
    orf_3p_tx <- as.integer(GenomicRanges::start(orf_3p_mapped))
    boundaries_exonic <-
      orf_5p_tx > 0L & orf_3p_tx > 0L &
      as.character(GenomicRanges::strand(orf_5p_mapped)) != "*" &
      as.character(GenomicRanges::strand(orf_3p_mapped)) != "*" &
      orf_5p_tx <= orf_3p_tx

    coding_transcript <- pair_tx_ids %in% names(cds_by_tx)
    coding_i <- which(coding_transcript)
    if (length(coding_i)) {
      coding_transcript[coding_i] <-
        S4Vectors::elementNROWS(cds_by_tx[pair_tx_ids[coding_i]]) > 0L
    }
    coding_i <- which(coding_transcript)

    cds_5p_tx <- rep(NA_integer_, length(qh))
    cds_3p_tx <- rep(NA_integer_, length(qh))

    if (length(coding_i)) {
      pair_cds <- cds_by_tx[pair_tx_ids[coding_i]]
      cds_5p_g <- ORFik::startSites(
        pair_cds, asGR = TRUE, keep.names = FALSE, is.sorted = FALSE
      )
      cds_3p_g <- ORFik::stopSites(
        pair_cds, asGR = TRUE, keep.names = FALSE, is.sorted = FALSE
      )
      cds_5p_mapped <- ORFik::pmapToTranscriptF(
        cds_5p_g, pair_tx[coding_i],
        x.is.sorted = TRUE, tx.is.sorted = FALSE,
        set.seqlengths = FALSE
      )
      cds_3p_mapped <- ORFik::pmapToTranscriptF(
        cds_3p_g, pair_tx[coding_i],
        x.is.sorted = TRUE, tx.is.sorted = FALSE,
        set.seqlengths = FALSE
      )
      cds_5p_tx[coding_i] <-
        as.integer(GenomicRanges::start(cds_5p_mapped))
      cds_3p_tx[coding_i] <-
        as.integer(GenomicRanges::start(cds_3p_mapped))
    }

    cds_projected <-
      !is.na(cds_5p_tx) & !is.na(cds_3p_tx) &
      cds_5p_tx > 0L & cds_3p_tx > 0L & cds_5p_tx <= cds_3p_tx
    classifiable <- boundaries_exonic & coding_transcript & cds_projected

    start_in_cds_frame <- rep(NA, length(qh))
    end_in_cds_frame <- rep(NA, length(qh))
    start_in_cds_frame[classifiable] <-
      (orf_5p_tx[classifiable] - cds_5p_tx[classifiable]) %% 3L == 0L
    end_in_cds_frame[classifiable] <-
      (orf_3p_tx[classifiable] - cds_3p_tx[classifiable]) %% 3L == 0L

    reference_orf_type <- rep(NA_character_, length(qh))
    reference_orf_type[boundaries_exonic & !coding_transcript] <-
      "varRNA-ORF"
    reference_orf_type[classifiable] <- .classify_orf_cds_geometry(
      orf_5p_tx[classifiable], orf_3p_tx[classifiable],
      cds_5p_tx[classifiable], cds_3p_tx[classifiable]
    )

    annotation_status <- dplyr::case_when(
      !boundaries_exonic ~ "ORF boundary not exonic on transcript",
      !coding_transcript ~ "compatible transcript has no annotated CDS",
      !cds_projected ~ "CDS could not be projected onto transcript",
      .default = "classified against annotated CDS"
    )

    pair_table <- tibble::tibble(
      orf_id = orf_ids[qh],
      transcript_id = pair_tx_ids,
      coding_transcript = coding_transcript,
      boundaries_exonic = boundaries_exonic,
      orf_5p_tx = orf_5p_tx,
      orf_3p_tx = orf_3p_tx,
      cds_5p_tx = cds_5p_tx,
      cds_3p_tx = cds_3p_tx,
      start_in_cds_frame = start_in_cds_frame,
      end_in_cds_frame = end_in_cds_frame,
      reference_orf_type = reference_orf_type,
      annotation_status = annotation_status
    )
  }

  valid_pairs <- pair_table |>
    dplyr::filter(boundaries_exonic, !is.na(reference_orf_type)) |>
    dplyr::mutate(.priority = match(reference_orf_type, type_priority))

  best_pair <- valid_pairs |>
    dplyr::arrange(orf_id, .priority, transcript_id) |>
    dplyr::group_by(orf_id) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      orf_id = orf_id,
      matched_transcript_id = transcript_id,
      reference_orf_type = reference_orf_type,
      reference_start_in_frame = start_in_cds_frame,
      reference_end_in_frame = end_in_cds_frame,
      reference_annotation_status = annotation_status
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

  orf_table <- tibble::tibble(
    orf_id = orf_ids,
    seqnames = as.character(GenomicRanges::seqnames(orfs)),
    start = GenomicRanges::start(orfs),
    end = GenomicRanges::end(orfs),
    strand = as.character(GenomicRanges::strand(orfs))
  ) |>
    dplyr::left_join(best_pair, by = "orf_id") |>
    dplyr::left_join(pair_summary, by = "orf_id") |>
    dplyr::mutate(
      reference_orf_type = dplyr::coalesce(
        reference_orf_type, "varRNA-ORF"
      ),
      reference_annotation_status = dplyr::coalesce(
        reference_annotation_status,
        "no transcript carries both ORF boundaries; fallback"
      ),
      n_compatible_transcripts = dplyr::coalesce(
        n_compatible_transcripts, 0L
      ),
      n_reference_orf_types = dplyr::coalesce(n_reference_orf_types, 0L),
      all_reference_orf_types = dplyr::coalesce(
        all_reference_orf_types, "varRNA-ORF"
      ),
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
    )

  list(orfs = orf_table, pairs = pair_table)
}
