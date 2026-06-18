# ==============================================================================
# orf_footprints.R
#
# Core data structure and constructor functions for ORF footprint objects.
#
# Three-level identity model (many-to-one upward):
#
#   transcript_id  →  chain_id  →  translon_id
#
# * translon_id  — stop + frame + strand + chrom. Biological translon.
#                  Grouping key only, NOT a primary key.
# * chain_id     — hash of: start + stop + strand + ordered internal intron
#                  boundaries between start and stop. Two calls sharing the
#                  same spliced path over the ORF span hash identically.
# * transcript_id — annotation label. Survives as set-valued metadata only.
#
# Output: GRangesList, one element per chain_id. mcols hold all per-chain
# metadata including list-column `compatible_transcripts`.
#
# Convention used throughout:
#   • Genomic coordinates are 1-based closed [start, end] (IRanges/GRanges).
#   • On the plus  strand: genomic start  <  genomic stop (start = A of start
#     codon, stop = last nt of stop codon).
#   • On the minus strand: genomic start  >  genomic stop, i.e. the A of the
#     start codon is at a HIGHER genomic coordinate than the stop codon.
#   • "Frame 0" = the first codon position is the very first nucleotide of the
#     ORF (i.e. the spliced 5'-most base on the coding strand).
#
# Dependencies: GenomicRanges, IRanges, S4Vectors, GenomicFeatures, dplyr,
#               tibble, digest, cli.
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports (used via :: throughout; listed here for NAMESPACE generation)
# ------------------------------------------------------------------------------
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @importFrom GenomicFeatures exonsBy
#' @importFrom dplyr tibble mutate filter group_by summarise n left_join
#'   bind_rows
#' @importFrom tibble tibble
#' @importFrom digest digest
#' @importFrom cli cli_abort cli_warn cli_inform
NULL


# ==============================================================================
# 1. check_txid_join_coverage
# ==============================================================================

#' Check Transcript ID Join Coverage Against a TxDb
#'
#' Run BEFORE any reconstruction step. Reports how many caller-supplied
#' `transcript_id`s resolve against the TxDb, with explicit handling of
#' GENCODE version suffixes (`ENST00000123456.7` → `ENST00000123456`).
#'
#' Match tiers (applied in order, mutually exclusive):
#' \enumerate{
#'   \item **exact** — id found verbatim in TxDb.
#'   \item **strip** — id found after removing `\\.\\d+$` suffix.
#'   \item **unresolvable** — not found at either tier.
#' }
#'
#' PAR-locus IDs (`_PAR_Y` suffix) and non-standard scaffold IDs
#' (seqnames not in the standard 1–22, X, Y, MT / chr1–22, chrX, chrY, chrM
#' set) are flagged separately within each tier.
#'
#' Nothing is dropped. The function is purely diagnostic.
#'
#' @param calls_df A data frame with at least a `transcript_id` character
#'   column. Additional columns are ignored.
#' @param txdb A TxDb object (e.g. from [GenomicFeatures::makeTxDbFromGFF()]).
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`summary`}{A [tibble::tibble()] with columns `tier`, `n`,
#'       `pct`, `n_par`, `n_nonstandard_scaffold`.}
#'     \item{`unresolvable_ids`}{Character vector of transcript_ids that could
#'       not be resolved at any tier.}
#'   }
#'
#' @export
check_txid_join_coverage <- function(calls_df, txdb) {
  if (!is.data.frame(calls_df)) {
    cli::cli_abort("{.arg calls_df} must be a data frame.")
  }
  if (!"transcript_id" %in% colnames(calls_df)) {
    cli::cli_abort("{.arg calls_df} must have a {.field transcript_id} column.")
  }
  if (!is(txdb, "TxDb")) {
    cli::cli_abort("{.arg txdb} must be a TxDb object.")
  }

  # Unique IDs from callers (work on unique set, report in caller-row terms)
  all_ids   <- calls_df$transcript_id
  unique_ids <- unique(all_ids)
  n_total    <- length(unique_ids)

  # All transcript names known to the TxDb.
  # AnnotationDbi::keys() with keytype="TXNAME" is the stable accessor.
  txdb_names <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
  txdb_names <- txdb_names[!is.na(txdb_names)]

  # Standard chromosome names (both UCSC and Ensembl styles accepted)
  .standard_chroms <- c(
    paste0("chr", c(1:22, "X", "Y", "M")),
    c(1:22, "X", "Y", "MT")
  )

  # Helper: is an id a PAR-locus entry?
  .is_par <- function(ids) grepl("_PAR_Y$", ids, ignore.case = FALSE)

  # Retrieve the seqname for a transcript id from the TxDb (vectorised).
  # AnnotationDbi::select() is the portable accessor across TxDb versions.
  .get_seqname <- function(ids, txdb) {
    res <- suppressMessages(AnnotationDbi::select(
      txdb,
      keys    = ids,
      columns = "TXCHROM",
      keytype = "TXNAME"
    ))
    seqnms <- res$TXCHROM
    names(seqnms) <- res$TXNAME
    seqnms[ids]
  }

  # ----- Tier 1: exact match --------------------------------------------------
  exact_mask <- unique_ids %in% txdb_names
  exact_ids  <- unique_ids[exact_mask]

  # ----- Tier 2: strip version suffix -----------------------------------------
  remaining       <- unique_ids[!exact_mask]
  stripped        <- sub("\\.\\d+$", "", remaining)
  strip_hit_mask  <- stripped %in% txdb_names
  strip_ids       <- remaining[strip_hit_mask]

  # ----- Tier 3: unresolvable -------------------------------------------------
  unresolvable_ids <- remaining[!strip_hit_mask]

  # ----- PAR / non-standard flags for resolved tiers -------------------------
  .flag_resolved <- function(ids, txdb, tier_label) {
    if (length(ids) == 0L) {
      return(tibble::tibble(
        tier = tier_label, n = 0L, pct = 0,
        n_par = 0L, n_nonstandard_scaffold = 0L
      ))
    }
    # For stripped IDs the TxDb name is the stripped version
    lookup <- if (tier_label == "strip") sub("\\.\\d+$", "", ids) else ids
    seqnms <- tryCatch(
      .get_seqname(lookup, txdb),
      error = function(e) rep(NA_character_, length(ids))
    )
    tibble::tibble(
      tier                  = tier_label,
      n                     = length(ids),
      pct                   = 100 * length(ids) / n_total,
      n_par                 = sum(.is_par(ids)),
      n_nonstandard_scaffold = sum(!is.na(seqnms) &
                                     !seqnms %in% .standard_chroms)
    )
  }

  summary_tbl <- dplyr::bind_rows(
    .flag_resolved(exact_ids,        txdb, "exact"),
    .flag_resolved(strip_ids,        txdb, "strip"),
    tibble::tibble(
      tier                   = "unresolvable",
      n                      = length(unresolvable_ids),
      pct                    = 100 * length(unresolvable_ids) / n_total,
      n_par                  = sum(.is_par(unresolvable_ids)),
      n_nonstandard_scaffold = NA_integer_
    )
  )

  if (length(unresolvable_ids) > 0) {
    cli::cli_warn(c(
      "!" = "{length(unresolvable_ids)} of {n_total} unique transcript_ids \\
             could not be resolved against the TxDb.",
      "i" = "Check {.code $unresolvable_ids} in the returned list."
    ))
  }

  list(
    summary          = summary_tbl,
    unresolvable_ids = unresolvable_ids
  )
}


# ==============================================================================
# 2. reconstruct_chain
# ==============================================================================

#' Reconstruct Spliced ORF Path from a Transcript
#'
#' Given genomic start, stop, strand, and a transcript ID, clips the
#' transcript's exons to the ORF span and returns only the exonic segments
#' between (and including) the start and stop positions.
#'
#' The returned GRanges is ALWAYS in 5′→3′ coding order:
#' * plus  strand: ascending genomic coordinate.
#' * minus strand: descending genomic coordinate (first element has the
#'   largest genomic position, i.e. the A of the start codon).
#'
#' @param start  Genomic coordinate of the first nt of the start codon (A of
#'   ATG). On minus strand this is numerically **greater** than `stop`.
#' @param stop   Genomic coordinate of the last nt of the stop codon. On minus
#'   strand this is numerically **less** than `start`.
#' @param strand Character `"+"` or `"-"`.
#' @param transcript_id Transcript name as it appears in `exons_by_tx` (after
#'   any version-stripping has been applied upstream).
#' @param exons_by_tx A named [GenomicRanges::GRangesList()] produced by
#'   `GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)`. Must be
#'   precomputed once and reused across calls.
#'
#' @return A [GenomicRanges::GRanges()] of clipped exon blocks in 5′→3′
#'   coding order, with a `exon_rank` integer metadata column (1 = first
#'   coding exon). Returns `NULL` (with a warning) if reconstruction fails.
#'
#' @export
reconstruct_chain <- function(start, stop, strand, transcript_id,
                               exons_by_tx) {
  # ------------------------------------------------------------------
  # Input guards
  # ------------------------------------------------------------------
  if (!strand %in% c("+", "-")) {
    cli::cli_abort("strand must be '+' or '-', got {.val {strand}}.")
  }
  if (!transcript_id %in% names(exons_by_tx)) {
    cli::cli_warn(
      "transcript_id {.val {transcript_id}} not found in exons_by_tx."
    )
    return(NULL)
  }

  exons <- exons_by_tx[[transcript_id]]

  # Ensure exons are on the expected strand
  # (TxDb exonsBy should already be strand-consistent, but guard anyway)
  if (!all(as.character(GenomicRanges::strand(exons)) == strand)) {
    cli::cli_warn(
      "Exons for {.val {transcript_id}} have mixed or unexpected strand. \\
       Skipping."
    )
    return(NULL)
  }

  # ------------------------------------------------------------------
  # Normalise: orf_left = min genomic coord, orf_right = max genomic coord
  # regardless of strand — we clip against the genomic interval.
  # ------------------------------------------------------------------
  orf_left  <- min(start, stop)   # numerically smaller genomic position
  orf_right <- max(start, stop)   # numerically larger  genomic position

  # ------------------------------------------------------------------
  # Select exons that overlap the ORF span [orf_left, orf_right].
  # An exon overlaps iff exon.start <= orf_right AND exon.end >= orf_left.
  # ------------------------------------------------------------------
  ex_start <- GenomicRanges::start(exons)
  ex_end   <- GenomicRanges::end(exons)

  overlapping <- which(ex_start <= orf_right & ex_end >= orf_left)

  if (length(overlapping) == 0L) {
    cli::cli_warn(
      "ORF [{orf_left}, {orf_right}] on {strand} does not overlap any exon \\
       of {.val {transcript_id}}. Caller/annotation mismatch?"
    )
    return(NULL)
  }

  chain <- exons[overlapping]

  # ------------------------------------------------------------------
  # Verify that start and stop each fall within an exon.
  # We check orf_left falls in the leftmost overlapping exon range, and
  # orf_right in the rightmost (after sorting by genomic position).
  # ------------------------------------------------------------------
  chain_sorted_asc <- chain[order(GenomicRanges::start(chain))]
  first_ex <- chain_sorted_asc[1]
  last_ex  <- chain_sorted_asc[length(chain_sorted_asc)]

  if (orf_left  < GenomicRanges::start(first_ex) ||
      orf_left  > GenomicRanges::end(first_ex)) {
    cli::cli_warn(
      "ORF left boundary {orf_left} is not within an exon of \\
       {.val {transcript_id}}. start/stop not within exon — \\
       caller/annotation mismatch."
    )
    return(NULL)
  }
  if (orf_right < GenomicRanges::start(last_ex) ||
      orf_right > GenomicRanges::end(last_ex)) {
    cli::cli_warn(
      "ORF right boundary {orf_right} is not within an exon of \\
       {.val {transcript_id}}. start/stop not within exon — \\
       caller/annotation mismatch."
    )
    return(NULL)
  }

  # ------------------------------------------------------------------
  # Clip first and last exon to the ORF boundaries.
  # ------------------------------------------------------------------
  chain_sorted_asc <- GenomicRanges::trim(chain_sorted_asc)

  # Clip the genomic-left exon: set its start to orf_left
  GenomicRanges::start(chain_sorted_asc[1]) <- orf_left

  # Clip the genomic-right exon: set its end to orf_right
  n_ex <- length(chain_sorted_asc)
  GenomicRanges::end(chain_sorted_asc[n_ex]) <- orf_right

  # ------------------------------------------------------------------
  # Order in 5′→3′ coding order:
  #   plus  strand → ascending  genomic coordinate (already sorted)
  #   minus strand → descending genomic coordinate
  # ------------------------------------------------------------------
  if (strand == "+") {
    # chain_sorted_asc is already in 5′→3′ order for plus strand
    coding_order <- seq_len(n_ex)
  } else {
    # minus strand: 5′ end is the highest genomic coordinate
    coding_order <- rev(seq_len(n_ex))
  }

  result <- chain_sorted_asc[coding_order]

  # Add exon rank (1 = first coding exon in 5′→3′ direction)
  GenomicRanges::mcols(result)$exon_rank <- seq_len(n_ex)

  result
}


# ==============================================================================
# 3. make_chain_id
# ==============================================================================

#' Compute a Deterministic Chain ID for an ORF Footprint
#'
#' A **chain_id** identifies the unique spliced path of an ORF — its codon
#' walk. Two calls with different `transcript_id`s that produce identical
#' exonic paths from start to stop will receive the **same** chain_id and
#' will be collapsed to a single footprint.
#'
#' ## What is hashed
#'
#' The canonical path representation is a character string of the form:
#' ```
#' strand:orf_left:orf_right:j1-j2:j3-j4:...
#' ```
#' where:
#' * `strand`    — `+` or `-`.
#' * `orf_left`  — numerically smaller genomic boundary (= start on `+`,
#'                 stop on `-`).
#' * `orf_right` — numerically larger  genomic boundary.
#' * `j1-j2`, `j3-j4`, … — **internal intron boundaries**, sorted
#'   numerically ascending, each formatted as `intron_start-intron_end`.
#'   An intron boundary is derived from consecutive exon blocks: if exon k
#'   ends at E and exon k+1 starts at S (both in ascending genomic order),
#'   the intron is `(E+1)-(S-1)`. Single-exon ORFs have no junction tokens.
#'
#' ## What is NOT hashed
#'
#' * `transcript_id` — by design. Cosmetic label differences must not
#'   create separate chains.
#' * UTR exons or exons entirely outside [orf_left, orf_right].
#' * Anything downstream of the stop codon.
#'
#' ## Identity guarantee
#'
#' Including both `orf_left`/`orf_right` (i.e. both start and stop) in the
#' hash ensures that transcripts sharing internal junctions but differing in
#' start codon position produce **different** chain_ids, reflecting genuinely
#' different codon walks at the N-terminus.
#'
#' @param footprint A [GenomicRanges::GRanges()] as returned by
#'   [reconstruct_chain()]: clipped exon blocks in any order (the function
#'   sorts internally). Must all share the same strand.
#'
#' @return A length-1 character string: a SHA-1 hex digest prefixed with
#'   `"chain_"`, e.g. `"chain_3a7f9b..."`.
#'
#' @export
make_chain_id <- function(footprint) {
  if (!is(footprint, "GRanges") || length(footprint) == 0L) {
    cli::cli_abort("{.arg footprint} must be a non-empty GRanges.")
  }

  strand_val <- unique(as.character(GenomicRanges::strand(footprint)))
  if (length(strand_val) != 1L || !strand_val %in% c("+", "-")) {
    cli::cli_abort(
      "All ranges in {.arg footprint} must share a single '+' or '-' strand."
    )
  }

  # Sort by ascending genomic position — strand-independent for boundary
  # extraction.  We always enumerate intron boundaries in ascending genomic
  # order so that the hash is strand-symmetric (same string regardless of
  # which end is "left").
  fp_asc <- footprint[order(GenomicRanges::start(footprint))]
  starts <- GenomicRanges::start(fp_asc)
  ends   <- GenomicRanges::end(fp_asc)

  # ORF genomic extent
  orf_left  <- starts[1]
  orf_right <- ends[length(ends)]

  # Internal intron boundaries: gaps between consecutive exon blocks.
  # Intron k runs from ends[k]+1 to starts[k+1]-1.
  # We record the boundary as "ends[k]+1 - starts[k+1]-1" (intron coords).
  n_ex <- length(fp_asc)
  if (n_ex > 1L) {
    intron_lefts  <- ends[-n_ex]  + 1L   # first base of each intron
    intron_rights <- starts[-1L]  - 1L   # last  base of each intron
    # Junctions already in ascending order since fp_asc is sorted
    junctions <- paste(intron_lefts, intron_rights, sep = "-")
    junction_str <- paste(junctions, collapse = ":")
  } else {
    # Single-exon ORF — no internal junctions
    junction_str <- ""
  }

  # Canonical string: strand:left:right:[junctions]
  canonical <- paste(strand_val, orf_left, orf_right, junction_str, sep = ":")

  # SHA-1 digest of the canonical string. SHA-1 is sufficient for collision
  # resistance at the scale of a single transcriptome (~1M unique chains).
  paste0("chain_", digest::digest(canonical, algo = "sha1", serialize = FALSE))
}


# ==============================================================================
# 4. compute_frame
# ==============================================================================

#' Compute the Reading Frame of an ORF Footprint
#'
#' Reading frame is defined relative to the chromosome coordinate system, as
#' used in GTF phase fields:
#'
#' **Frame 0** — the first nucleotide of the ORF (5′ end on the coding strand)
#' is at genomic position `pos` where `(pos - 1) %% 3 == 0` relative to the
#' chromosome origin. More practically: the start codon occupies the first
#' three spliced nucleotides (positions 1–3 in 1-based spliced coordinates).
#'
#' Frame is computed from the genomic position of the first coding base:
#' * plus  strand: `(start_pos - 1) %% 3`
#' * minus strand: `(end_pos   - 1) %% 3`  where `end_pos` is the last
#'   genomic nt of the start codon (numerically largest on the minus strand).
#'
#' If a `phase` column is present in `mcols(footprint)` (from a TxDb-derived
#' CDS GRanges) the GTF phase of the first CDS block is used directly:
#' `phase` ∈ {0, 1, 2} where 0 means in-frame.
#'
#' @param footprint A [GenomicRanges::GRanges()] of clipped exon blocks as
#'   returned by [reconstruct_chain()].
#'
#' @return Integer 0, 1, or 2.
#'
#' @export
compute_frame <- function(footprint) {
  if (!is(footprint, "GRanges") || length(footprint) == 0L) {
    cli::cli_abort("{.arg footprint} must be a non-empty GRanges.")
  }

  strand_val <- unique(as.character(GenomicRanges::strand(footprint)))
  if (length(strand_val) != 1L) {
    cli::cli_abort("All blocks must share a single strand.")
  }

  # Prefer explicit phase annotation if present (e.g. from TxDb CDS ranges)
  if ("phase" %in% colnames(GenomicRanges::mcols(footprint))) {
    # Sort to get the 5′-most block, then read its phase
    if (strand_val == "+") {
      first_block <- footprint[which.min(GenomicRanges::start(footprint))]
    } else {
      # minus strand: 5′ end = highest genomic coord
      first_block <- footprint[which.max(GenomicRanges::end(footprint))]
    }
    phase_val <- GenomicRanges::mcols(first_block)$phase
    if (!is.na(phase_val) && phase_val %in% c(0L, 1L, 2L)) {
      return(as.integer(phase_val))
    }
    # Fall through if phase is NA or out of range
  }

  # Compute from the genomic position of the first coding nucleotide.
  # Convention: frame = (pos - 1) %% 3, where pos is 1-based genomic.
  if (strand_val == "+") {
    # plus strand: first coding nt = start of the leftmost exon block
    first_pos <- min(GenomicRanges::start(footprint))
  } else {
    # minus strand: first coding nt = end (highest coord) of rightmost block
    # because reading proceeds right→left on the minus strand
    first_pos <- max(GenomicRanges::end(footprint))
  }

  as.integer((first_pos - 1L) %% 3L)
}


# ==============================================================================
# 5. make_translon_id / make_start_id
# ==============================================================================

#' Make a Translon ID
#'
#' The translon is defined by stop codon position + reading frame + strand +
#' chromosome. It is a **grouping key**, not a primary key — multiple chains
#' (different start codons, different splicing upstream of a shared stop) can
#' share a translon_id.
#'
#' Format: `chrom:strand:stop_pos:fN`
#' e.g.    `chr1:+:1234567:f0`
#'
#' `stop_pos` is always the **last genomic nucleotide of the stop codon**:
#' * plus  strand: the numerically largest position in the footprint.
#' * minus strand: the numerically smallest position.
#'
#' @param footprint A [GenomicRanges::GRanges()] as from [reconstruct_chain()].
#' @return Length-1 character string.
#' @export
make_translon_id <- function(footprint) {
  .check_footprint(footprint)

  strand_val <- unique(as.character(GenomicRanges::strand(footprint)))
  chrom      <- unique(as.character(GenomicRanges::seqnames(footprint)))

  # Stop codon: 3′ end of the ORF on the coding strand.
  # plus  strand → rightmost genomic position
  # minus strand → leftmost  genomic position
  stop_pos <- if (strand_val == "+") {
    max(GenomicRanges::end(footprint))
  } else {
    min(GenomicRanges::start(footprint))
  }

  frame <- compute_frame(footprint)

  paste0(chrom, ":", strand_val, ":", stop_pos, ":f", frame)
}


#' Make a Start ID
#'
#' Extends [make_translon_id()] by additionally encoding the start codon
#' position. Useful for grouping footprints that share a start but differ in
#' their internal path (e.g. alternative splicing between start and stop).
#'
#' Format: `chrom:strand:start_pos:stop_pos:fN`
#'
#' `start_pos` is the A of the start codon:
#' * plus  strand: numerically smallest genomic position in the footprint.
#' * minus strand: numerically largest.
#'
#' @param footprint A [GenomicRanges::GRanges()] as from [reconstruct_chain()].
#' @return Length-1 character string.
#' @export
make_start_id <- function(footprint) {
  .check_footprint(footprint)

  strand_val <- unique(as.character(GenomicRanges::strand(footprint)))
  chrom      <- unique(as.character(GenomicRanges::seqnames(footprint)))

  start_pos <- if (strand_val == "+") {
    min(GenomicRanges::start(footprint))
  } else {
    max(GenomicRanges::end(footprint))
  }

  stop_pos <- if (strand_val == "+") {
    max(GenomicRanges::end(footprint))
  } else {
    min(GenomicRanges::start(footprint))
  }

  frame <- compute_frame(footprint)

  paste0(chrom, ":", strand_val, ":", start_pos, ":", stop_pos, ":f", frame)
}

# Internal guard shared by make_translon_id and make_start_id
.check_footprint <- function(fp) {
  if (!is(fp, "GRanges") || length(fp) == 0L)
    cli::cli_abort("{.arg footprint} must be a non-empty GRanges.")
  if (length(unique(as.character(GenomicRanges::seqnames(fp)))) != 1L)
    cli::cli_abort("All blocks must be on the same chromosome.")
  if (length(unique(as.character(GenomicRanges::strand(fp)))) != 1L)
    cli::cli_abort("All blocks must share a single strand.")
}


# ==============================================================================
# 6. parse_orf_footprints
# ==============================================================================

#' Parse ORF Calls into a Chain-Level GRangesList of Footprints
#'
#' Top-level constructor. Accepts a data frame of ORF calls (potentially from
#' multiple callers), reconstructs the spliced path for each call, groups by
#' `chain_id`, and returns one [GenomicRanges::GRangesList()] element per
#' distinct chain.
#'
#' Collapse semantics:
#' * Two calls producing the same spliced ORF path → **one** footprint.
#' * `compatible_transcripts` in `mcols`: all transcript IDs whose paths
#'   collapsed here (set-valued; stored as a `CharacterList`).
#' * `called_by` in `mcols`: named list mapping caller → transcript_id it
#'   used for this chain.
#' * Calls whose reconstruction failed are retained with
#'   `chain_reconstructed = FALSE` and a synthetic GRanges placeholder
#'   (the genomic interval start–stop with no exon structure).
#'
#' @param calls_df A data frame with columns:
#'   \describe{
#'     \item{`start`}{Integer. Genomic start (A of start codon; plus strand:
#'       smaller coord; minus strand: larger coord).}
#'     \item{`stop`}{Integer. Genomic stop (last nt of stop codon).}
#'     \item{`strand`}{Character `"+"` or `"-"`.}
#'     \item{`transcript_id`}{Character. GENCODE transcript ID (with or without
#'       version suffix). Will be matched against `exons_by_tx` names; version
#'       stripping is attempted automatically.}
#'     \item{`caller`}{Character. Name of the calling tool
#'       (e.g. `"RiboTIE"`, `"RiboCode"`, `"ORFquant"`).}
#'     \item{`run`}{Character, **optional**. Sample/run identifier (e.g.
#'       `"RT_HEK293_3"`). When present, aggregated into a `called_by_run`
#'       list-column used by [filter_chains_by_junctions()]. Rows with
#'       `NA` run (e.g. from all-transcript enumeration) are silently excluded
#'       from `called_by_run`.}
#'   }
#' @param txdb A TxDb object.
#'
#' @return A named [GenomicRanges::GRangesList()] where names are `chain_id`
#'   values, and `mcols` contains:
#'   \describe{
#'     \item{`chain_id`}{Character.}
#'     \item{`translon_id`}{Character.}
#'     \item{`start_id`}{Character.}
#'     \item{`n_exons`}{Integer. Number of exon blocks in the chain.}
#'     \item{`single_exon`}{Logical.}
#'     \item{`frame`}{Integer 0/1/2.}
#'     \item{`compatible_transcripts`}{[IRanges::CharacterList()].}
#'     \item{`called_by`}{[S4Vectors::SimpleList()] of named character vectors,
#'       caller → transcript_ids.}
#'     \item{`called_by_run`}{[IRanges::CharacterList()]. Run names that called
#'       this chain. Only present when `calls_df` has a `run` column.}
#'     \item{`n_calls`}{Integer. Number of raw calls collapsed here.}
#'     \item{`chain_reconstructed`}{Logical. FALSE if reconstruction failed.}
#'   }
#'
#' @export
parse_orf_footprints <- function(calls_df, txdb) {
  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------
  required_cols <- c("start", "stop", "strand", "transcript_id", "caller")
  missing_cols  <- setdiff(required_cols, colnames(calls_df))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      "{.arg calls_df} is missing required columns: \\
       {.field {missing_cols}}."
    )
  }
  if (!is(txdb, "TxDb")) {
    cli::cli_abort("{.arg txdb} must be a TxDb object.")
  }

  n_calls <- nrow(calls_df)
  cli::cli_inform(c("i" = "Processing {n_calls} ORF call{?s}..."))

  # ------------------------------------------------------------------
  # Precompute exonsBy ONCE — never inside any per-row operation.
  # ------------------------------------------------------------------
  cli::cli_inform(c("i" = "Precomputing exonsBy(txdb) ..."))
  exons_by_tx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)

  # ------------------------------------------------------------------
  # Attempt version-suffix resolution for transcript IDs not found verbatim.
  # ------------------------------------------------------------------
  txdb_names   <- names(exons_by_tx)
  raw_tx_ids   <- calls_df$transcript_id
  needs_strip  <- !raw_tx_ids %in% txdb_names
  stripped_ids <- sub("\\.\\d+$", "", raw_tx_ids)
  resolved_ids <- ifelse(
    needs_strip & stripped_ids %in% txdb_names,
    stripped_ids,
    raw_tx_ids
  )
  calls_df$resolved_tx_id <- resolved_ids

  # ------------------------------------------------------------------
  # Reconstruct chain for each call (vectorised over rows via lapply)
  # ------------------------------------------------------------------
  cli::cli_inform(c("i" = "Reconstructing spliced paths..."))

  reconstructions <- lapply(seq_len(n_calls), function(i) {
    tryCatch(
      reconstruct_chain(
        start         = calls_df$start[i],
        stop          = calls_df$stop[i],
        strand        = calls_df$strand[i],
        transcript_id = calls_df$resolved_tx_id[i],
        exons_by_tx   = exons_by_tx
      ),
      error = function(e) NULL
    )
  })

  chain_reconstructed <- !vapply(reconstructions, is.null, logical(1))

  # ------------------------------------------------------------------
  # Compute chain_id for successfully reconstructed calls;
  # use a placeholder key for failures.
  # ------------------------------------------------------------------
  chain_ids <- character(n_calls)

  for (i in seq_len(n_calls)) {
    if (chain_reconstructed[i]) {
      chain_ids[i] <- make_chain_id(reconstructions[[i]])
    } else {
      # Deterministic fallback key encoding the raw coordinates + strand
      # (NOT the transcript_id — preserves the identity semantics even for
      # failed reconstructions).
      fallback_str <- paste(
        calls_df$strand[i],
        min(calls_df$start[i], calls_df$stop[i]),
        max(calls_df$start[i], calls_df$stop[i]),
        sep = ":"
      )
      chain_ids[i] <- paste0(
        "chain_failed_",
        digest::digest(fallback_str, algo = "sha1", serialize = FALSE)
      )
    }
  }

  calls_df$chain_id           <- chain_ids
  calls_df$chain_reconstructed <- chain_reconstructed

  # ------------------------------------------------------------------
  # Group by chain_id and aggregate
  # ------------------------------------------------------------------
  cli::cli_inform(c("i" = "Collapsing {n_calls} calls to unique chains..."))

  unique_chains <- unique(chain_ids)
  n_chains      <- length(unique_chains)

  cli::cli_inform(c(
    "v" = "Found {n_chains} unique chain{?s} from {n_calls} call{?s}."
  ))

  # Build the GRangesList: one element per chain
  chain_grl    <- vector("list", n_chains)
  chain_mcols  <- vector("list", n_chains)

  for (k in seq_len(n_chains)) {
    cid    <- unique_chains[k]
    idx    <- which(chain_ids == cid)

    # Representative call index: first successfully reconstructed call in this
    # group, falling back to any call if all failed.
    rep_idx <- idx[chain_reconstructed[idx]][1]
    if (is.na(rep_idx)) rep_idx <- idx[1]

    rep_fp <- reconstructions[[rep_idx]]
    ok     <- chain_reconstructed[rep_idx]

    # Aggregate across all calls in this group
    comp_tx  <- unique(calls_df$resolved_tx_id[idx])
    n_in_grp <- length(idx)

    # called_by: named list of caller → transcript_id(s) it used for this chain
    caller_vec   <- calls_df$caller[idx]
    tx_vec       <- calls_df$resolved_tx_id[idx]
    called_by_df <- split(tx_vec, caller_vec)

    # called_by_run: unique run names that produced this chain (NA rows excluded)
    # Only populated when calls_df has a run column.
    if ("run" %in% colnames(calls_df)) {
      run_vec      <- calls_df$run[idx]
      called_by_run_vec <- unique(run_vec[!is.na(run_vec)])
    } else {
      called_by_run_vec <- character(0L)
    }

    if (ok) {
      chain_grl[[k]]   <- rep_fp
      n_ex             <- length(rep_fp)
      frame_val        <- compute_frame(rep_fp)
      translon_val     <- make_translon_id(rep_fp)
      start_id_val     <- make_start_id(rep_fp)
    } else {
      # Placeholder: raw genomic interval as a single-block GRanges
      s_left  <- min(calls_df$start[rep_idx], calls_df$stop[rep_idx])
      s_right <- max(calls_df$start[rep_idx], calls_df$stop[rep_idx])
      rep_chrom <- NA_character_  # seqname not reliably known from calls_df
      rep_fp_placeholder <- GenomicRanges::GRanges(
        seqnames = if ("chrom" %in% colnames(calls_df))
          calls_df$chrom[rep_idx]
        else
          "unknown",
        ranges   = IRanges::IRanges(start = s_left, end = s_right),
        strand   = calls_df$strand[rep_idx]
      )
      chain_grl[[k]] <- rep_fp_placeholder
      n_ex           <- NA_integer_
      frame_val      <- NA_integer_
      translon_val   <- NA_character_
      start_id_val   <- NA_character_
    }

    chain_mcols[[k]] <- list(
      chain_id               = cid,
      translon_id            = translon_val,
      start_id               = start_id_val,
      n_exons                = as.integer(n_ex),
      single_exon            = if (!is.na(n_ex)) n_ex == 1L else NA,
      frame                  = as.integer(frame_val),
      compatible_transcripts = comp_tx,
      called_by              = called_by_df,
      called_by_run          = called_by_run_vec,
      n_calls                = as.integer(n_in_grp),
      chain_reconstructed    = ok
    )
  }

  # ------------------------------------------------------------------
  # Assemble GRangesList with mcols
  # ------------------------------------------------------------------
  grl <- GenomicRanges::GRangesList(chain_grl)
  names(grl) <- unique_chains

  # Build mcols DataFrame — list-columns need explicit construction
  has_run_col <- "run" %in% colnames(calls_df)

  mc_df <- S4Vectors::DataFrame(
    chain_id            = vapply(chain_mcols, `[[`, character(1), "chain_id"),
    translon_id         = vapply(chain_mcols, `[[`, character(1), "translon_id"),
    start_id            = vapply(chain_mcols, `[[`, character(1), "start_id"),
    n_exons             = vapply(chain_mcols, `[[`, integer(1),   "n_exons"),
    single_exon         = vapply(chain_mcols, `[[`, logical(1),   "single_exon"),
    frame               = vapply(chain_mcols, `[[`, integer(1),   "frame"),
    n_calls             = vapply(chain_mcols, `[[`, integer(1),   "n_calls"),
    chain_reconstructed = vapply(chain_mcols, `[[`, logical(1),   "chain_reconstructed"),
    compatible_transcripts = IRanges::CharacterList(
      lapply(chain_mcols, `[[`, "compatible_transcripts")
    ),
    called_by = S4Vectors::SimpleList(
      lapply(chain_mcols, `[[`, "called_by")
    )
  )

  if (has_run_col) {
    mc_df$called_by_run <- IRanges::CharacterList(
      lapply(chain_mcols, `[[`, "called_by_run")
    )
  }

  S4Vectors::mcols(grl) <- mc_df

  grl
}


# ==============================================================================
# 7. filter_chains_by_junctions
# ==============================================================================

#' Filter Chain-Level Footprints by BAM Junction Evidence
#'
#' Retains only chains whose internal intron boundaries are supported by
#' junction-spanning reads in one or more BAM files. Single-exon chains and
#' chains with failed reconstruction pass unconditionally (they have no
#' junctions to check).
#'
#' @section Algorithm:
#' 1. Extract all internal intron GRanges from multi-exon chains.
#' 2. Read junction-spanning reads from each BAM restricted to those regions
#'    via [GenomicAlignments::summarizeJunctions()].
#' 3. Match observed junctions to chain junctions by exact position
#'    ([GenomicRanges::findOverlaps()] with `type = "equal"`).
#' 4. A chain passes if **every** one of its junctions has ≥ `min_reads`
#'    support in the relevant BAMs.
#'
#' @param footprints A named [GenomicRanges::GRangesList()] as returned by
#'   [parse_orf_footprints()].
#' @param bam_files A **named** character vector of BAM file paths. Names must
#'   match run identifiers in `mcols(footprints)$called_by_run` when
#'   `require_same_sample = TRUE`.
#' @param min_reads Integer. Minimum number of junction-spanning reads required
#'   per junction. Default `2L`.
#' @param require_same_sample Logical. If `TRUE`, a chain is assessed only
#'   against BAMs from the runs that called it (via `called_by_run` in
#'   `mcols`). Requires `parse_orf_footprints()` to have been run with a
#'   `calls_df` that had a `run` column. If no matching BAMs are found for a
#'   chain, it is kept (cannot be assessed). Default `FALSE` (all BAMs pooled).
#'
#' @return A filtered [GenomicRanges::GRangesList()] with the same structure as
#'   the input. A `junction_reads` integer list-column is added to `mcols`,
#'   recording the per-junction read counts used for the decision
#'   (`NA` for single-exon or failed-reconstruction chains).
#'
#' @importFrom GenomicAlignments readGAlignments summarizeJunctions
#' @importFrom Rsamtools ScanBamParam
#' @export
filter_chains_by_junctions <- function(footprints,
                                        bam_files,
                                        min_reads           = 2L,
                                        require_same_sample = FALSE) {

  if (!is(footprints, "GRangesList") || length(footprints) == 0L) {
    cli::cli_abort("{.arg footprints} must be a non-empty GRangesList.")
  }
  if (!is.character(bam_files) || is.null(names(bam_files))) {
    cli::cli_abort("{.arg bam_files} must be a named character vector.")
  }
  missing_bams <- bam_files[!file.exists(bam_files)]
  if (length(missing_bams) > 0L) {
    cli::cli_abort("BAM file{?s} not found: {.path {missing_bams}}")
  }
  if (require_same_sample &&
      !"called_by_run" %in% colnames(S4Vectors::mcols(footprints))) {
    cli::cli_abort(
      "{.arg require_same_sample} requires a {.field called_by_run} column in \\
       {.code mcols(footprints)}. Re-run {.fn parse_orf_footprints} with a \\
       {.field run} column in {.arg calls_df}."
    )
  }

  is_single <- S4Vectors::mcols(footprints)$single_exon
  is_failed <- !S4Vectors::mcols(footprints)$chain_reconstructed

  # Indices of chains that need junction checking
  multi_idx <- which(!is_single & !is_failed)

  # Initialise junction_reads list-column (NA = not assessed)
  junction_reads_list <- vector("list", length(footprints))

  if (length(multi_idx) == 0L) {
    cli::cli_inform("All reconstructed chains are single-exon; nothing to filter.")
    S4Vectors::mcols(footprints)$junction_reads <-
      S4Vectors::SimpleList(junction_reads_list)
    return(footprints)
  }

  # ------------------------------------------------------------------
  # Extract internal intron GRanges for all multi-exon chains
  # ------------------------------------------------------------------
  junction_list <- lapply(multi_idx, function(i) {
    fp    <- footprints[[i]]
    cid   <- S4Vectors::mcols(footprints)$chain_id[i]

    # Sort ascending for intron coordinate extraction (strand-independent)
    fp_asc <- fp[order(GenomicRanges::start(fp))]
    n_ex   <- length(fp_asc)

    intron_starts <- GenomicRanges::end(fp_asc)[-n_ex]        + 1L
    intron_ends   <- GenomicRanges::start(fp_asc)[-1L]        - 1L

    gr <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(fp_asc)[1],
      ranges   = IRanges::IRanges(intron_starts, intron_ends),
      strand   = GenomicRanges::strand(fp_asc)[1]
    )
    gr$chain_id   <- cid
    gr$chain_idx  <- i          # back-reference to footprints index
    gr$junction_n <- seq_len(n_ex - 1L)
    gr
  })

  all_junctions <- do.call(c, junction_list)

  # ------------------------------------------------------------------
  # Accumulate junction read counts from each BAM
  # Matrix: rows = all_junctions, cols = BAM runs
  # Reads are restricted to junction regions for efficiency.
  # ------------------------------------------------------------------
  n_junc <- length(all_junctions)
  counts_mat <- matrix(
    0L,
    nrow     = n_junc,
    ncol     = length(bam_files),
    dimnames = list(NULL, names(bam_files))
  )

  # Query regions: the intron intervals themselves — junction-spanning reads
  # have CIGAR N operations that span these coordinates, so their full
  # alignment range overlaps the intron region.
  query_regions <- GenomicRanges::reduce(all_junctions, ignore.strand = TRUE)
  bam_param     <- Rsamtools::ScanBamParam(which = query_regions)

  cli::cli_inform(c(
    "i" = "Querying {length(bam_files)} BAM file{?s} for junction evidence \\
           across {n_junc} junction{?s}..."
  ))

  for (b in seq_along(bam_files)) {
    ga <- GenomicAlignments::readGAlignments(
      bam_files[b],
      param = bam_param
    )
    if (length(ga) == 0L) next

    obs_junc <- GenomicAlignments::summarizeJunctions(
      ga,
      ignore.strand = FALSE
    )
    if (length(obs_junc) == 0L) next

    # Exact-position match: intron start and end must be identical
    hits <- GenomicRanges::findOverlaps(
      all_junctions, obs_junc,
      type         = "equal",
      ignore.strand = FALSE
    )
    counts_mat[S4Vectors::queryHits(hits), b] <-
      GenomicRanges::score(obs_junc)[S4Vectors::subjectHits(hits)]
  }

  # ------------------------------------------------------------------
  # Decide pass/fail per chain; record per-junction counts
  # ------------------------------------------------------------------
  chain_passes <- logical(length(multi_idx))

  for (ki in seq_along(multi_idx)) {
    i   <- multi_idx[ki]
    cid <- S4Vectors::mcols(footprints)$chain_id[i]

    j_rows <- which(all_junctions$chain_id == cid)   # rows in counts_mat

    if (require_same_sample) {
      run_names <- as.character(
        S4Vectors::mcols(footprints)$called_by_run[[i]]
      )
      bam_cols  <- intersect(run_names, colnames(counts_mat))

      if (length(bam_cols) == 0L) {
        # No matching BAMs — cannot assess; keep the chain
        chain_passes[ki]        <- TRUE
        junction_reads_list[[i]] <- rep(NA_integer_, length(j_rows))
        next
      }
      per_junc <- rowSums(counts_mat[j_rows, bam_cols, drop = FALSE])
    } else {
      per_junc <- rowSums(counts_mat[j_rows, , drop = FALSE])
    }

    junction_reads_list[[i]] <- as.integer(per_junc)
    chain_passes[ki]         <- all(per_junc >= min_reads)
  }

  # ------------------------------------------------------------------
  # Apply filter and report
  # ------------------------------------------------------------------
  keep          <- rep(TRUE, length(footprints))
  keep[multi_idx] <- chain_passes

  n_removed <- sum(!keep)
  cli::cli_inform(c(
    "v" = "Junction filter: removed {n_removed} of {length(multi_idx)} \\
           multi-exon chain{?s}.",
    "i" = "Retained {sum(keep)} chain{?s} \\
           ({sum(is_single, na.rm = TRUE)} single-exon, \\
           {sum(keep & !is_single, na.rm = TRUE)} multi-exon)."
  ))

  result <- footprints[keep]
  S4Vectors::mcols(result)$junction_reads <-
    S4Vectors::SimpleList(junction_reads_list[keep])

  result
}
