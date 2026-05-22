#' Check Sequence Level Compatibility Between Genomic Objects
#'
#' Validates that two genomic objects (e.g., TxDb and BSgenome) have compatible
#' sequence levels and naming styles. This prevents mismatches between chromosome
#' naming conventions (e.g., "chr1" vs "1").
#'
#' @param x,y Genomic objects with sequence levels (e.g., TxDb, BSgenome, GRanges).
#'   Must support \\code{seqlevelsStyle()} and \\code{seqlevels()} methods from
#'   GenomeInfoDb.
#'
#' @return NULL (invisibly). Function is called for side effects (error on mismatch).
#'
#' @details
#' The function checks:\n#' \\enumerate{\n#'   \\item \\strong{Sequence style}: Ensures both objects use the same naming\n#'     convention (UCSC: \"chr1\", Ensembl: \"1\", NCBI: \"NC_000001.11\", etc.)\n#'   \\item \\strong{Sequence names}: Verifies that sequence names are compatible\n#'     (though exact match is not required)\n#' }\n#'
#' Common mismatches occur when mixing:\n#' \\itemize{\n#'   \\item UCSC genome (BSgenome.Hsapiens.UCSC.hg38) with Ensembl GTF\n#'   \\item Chromosome-only GTF with whole-genome alignment including contigs\n#' }\n#'
#' @keywords internal
#' @export
#'
#' @examples
#' \\dontrun{\n#' library(BSgenome.Hsapiens.UCSC.hg38)\n#' library(GenomicFeatures)\n#' \n#' # Load GTF and genome\n#' txdb <- makeTxDbFromGFF(\"gencode.gtf\")\n#' genome <- BSgenome.Hsapiens.UCSC.hg38\n#' \n#' # Check compatibility (will error if incompatible)\n#' check_seq_levels(txdb, genome)\n#' }\n#'\n#' @seealso \\code{GenomeInfoDb::seqlevelsStyle} for viewing/converting styles
check_seq_levels <- function(x, y){

  if(!seqlevelsStyle(x) == seqlevelsStyle(y)){
    cli::cli_abort(
      "{.arg {arg}} must have the same seq levels style, see GenomeInfoDb::seqlevelsStyle()"
    )
  }

  if (setequal(seqlevels(x), seqlevels(y))){
    cli::cli_abort(
      "{.arg {arg}} must have the same seq levels"
    )
  }

}

#' Assign Strand-Aware Exon Ranks to Genomic Ranges
#'
#' Adds exon rank metadata to each range in a GRangesList, with directionality
#' based on strand. For positive strand, ranks 5' to 3' (1, 2, 3...). For
#' negative strand, ranks 3' to 5' (reverse order) so exon rank 1 is always
#' the 5'-most exon in transcript orientation.
#'
#' @param grl GRangesList where each element represents a transcript or ORF
#'   with multiple exons. Strand information must be present.
#'
#' @return A GRangesList with the same structure as input, but with an additional
#'   \\code{exon_rank} column in the metadata. Ranges are sorted by exon_rank
#'   within each element.
#'
#' @details
#' The function processes each element of the GRangesList independently:
#' \\itemize{
#'   \\item \\strong{Positive strand (+)}: Ranks 1, 2, 3... following genomic
#'     coordinates (left to right)
#'   \\item \\strong{Negative strand (-)}: Ranks from highest to lowest coordinate
#'     so rank 1 is rightmost (which is 5' in transcript orientation)
#'   \\item \\strong{Unstranded (*)}: Treated as positive strand
#' }
#'
#' This ensures exon rank 1 always refers to the first exon encountered during
#' transcription/translation, regardless of strand.
#'
#' @export
#'
#' @examples
#' library(GenomicRanges)
#' 
#' # Positive strand transcript
#' tx_pos <- GRanges(\"chr1\", 
#'                   IRanges(c(100, 200, 300), c(150, 250, 350)),
#'                   strand = \"+\")
#' grl_pos <- GRangesList(gene1 = tx_pos)
#' ranked_pos <- rank_exons(grl_pos)
#' ranked_pos$gene1$exon_rank  # 1, 2, 3
#' 
#' # Negative strand transcript
#' tx_neg <- GRanges(\"chr1\",
#'                   IRanges(c(100, 200, 300), c(150, 250, 350)),
#'                   strand = \"-\")
#' grl_neg <- GRangesList(gene1 = tx_neg)
#' ranked_neg <- rank_exons(grl_neg)
#' ranked_neg$gene1$exon_rank  # 3, 2, 1 (then sorted to 1, 2, 3 by position)
#' start(ranked_neg$gene1)     # 300, 200, 100 (sorted by rank)
#'
#' @seealso \\code{\\link{annotate_orf_isoforms}} which uses this function internally
rank_exons <- function(grl){
  # endoapply preserves GRangesList structure (length, names, type) when
  # applying a function to each element — safer than lapply + GRangesList().
  endoapply(grl, function(x){
    if (length(x) == 0) return(x)
    str <- as.character(strand(x))
    if (all(str == "-")) {
      x$exon_rank <- rev(seq_along(x))
    } else {
      x$exon_rank <- seq_along(x)
    }
    x[order(x$exon_rank, decreasing = FALSE)]
  })
}

#' Get downstream nucleotides after an ORF within a transcript
#'
#' Retrieves the \code{n} nucleotides immediately 3' of an ORF within the
#' spliced structure of its transcript. Used to check whether the next codon
#' after the ORF is a stop codon. Correctly handles exon-intron boundaries:
#' if the ORF ends in one exon and the downstream nucleotides span into the
#' next exon, both segments are returned as a multi-range GRanges.
#'
#' @param orf_exons GRanges of the ORF exonic segments (already intersected
#'   with the transcript exons).
#' @param tx_exons GRanges of the parent transcript's exons.
#' @param n integer, number of nucleotides to retrieve. Default 3 (one codon).
#'
#' @return A GRanges with the genomic positions of the \code{n} downstream
#'   nucleotides, or an empty GRanges if the transcript ends before \code{n}
#'   nucleotides are available.
#'
#' @keywords internal
.get_downstream_nt <- function(orf_exons, tx_exons, n = 3L) {
  # Strip metadata so exons carried through the else-branch below don't inherit
  # exon_rank (or any other mcols) from annotations$transcripts.  Without this,
  # mixing bare GRanges() results (no mcols) with subsets of tx_exons (may have
  # exon_rank, possibly NA) causes do.call(c, result) to produce a GRanges with
  # an exon_rank column containing NAs, which triggers ORFik's check.
  mcols(orf_exons) <- NULL
  mcols(tx_exons)  <- NULL
  strnd <- as.character(strand(orf_exons)[1L])

  if (strnd %in% c("+", "*")) {
    orf_3prime <- max(end(orf_exons))
    # Transcript exons that extend beyond the ORF 3' end
    ds_exons <- tx_exons[end(tx_exons) > orf_3prime]
    if (length(ds_exons) == 0L) return(GRanges())
    ds_exons <- sort(ds_exons)
    # Clip any exon that overlaps the ORF to start just after the ORF ends
    start(ds_exons) <- pmax(start(ds_exons), orf_3prime + 1L)
    ds_exons <- ds_exons[width(ds_exons) > 0L]
    if (length(ds_exons) == 0L) return(GRanges())

    remaining <- n
    result    <- list()
    for (i in seq_along(ds_exons)) {
      w <- width(ds_exons[i])
      if (w >= remaining) {
        result[[length(result) + 1L]] <- GRanges(
          seqnames(ds_exons[i]),
          IRanges(start(ds_exons[i]), start(ds_exons[i]) + remaining - 1L),
          strand = strand(ds_exons[i])
        )
        remaining <- 0L
        break
      } else {
        result[[length(result) + 1L]] <- ds_exons[i]
        remaining <- remaining - w
      }
    }
    if (remaining > 0L) return(GRanges())   # transcript ends before n nt
    do.call(c, result)

  } else {   # "-" strand: 3' end of ORF is at the lowest genomic coordinate
    orf_3prime <- min(start(orf_exons))
    ds_exons   <- tx_exons[start(tx_exons) < orf_3prime]
    if (length(ds_exons) == 0L) return(GRanges())
    end(ds_exons) <- pmin(end(ds_exons), orf_3prime - 1L)
    ds_exons <- ds_exons[width(ds_exons) > 0L]
    if (length(ds_exons) == 0L) return(GRanges())
    # Reverse: highest genomic position first = 5'→3' in transcript orientation
    ds_exons <- rev(sort(ds_exons))

    remaining <- n
    result    <- list()
    for (i in seq_along(ds_exons)) {
      w <- width(ds_exons[i])
      if (w >= remaining) {
        result[[length(result) + 1L]] <- GRanges(
          seqnames(ds_exons[i]),
          IRanges(end(ds_exons[i]) - remaining + 1L, end(ds_exons[i])),
          strand = strand(ds_exons[i])
        )
        remaining <- 0L
        break
      } else {
        result[[length(result) + 1L]] <- ds_exons[i]
        remaining <- remaining - w
      }
    }
    if (remaining > 0L) return(GRanges())
    do.call(c, result)
  }
}

# Returns the 0-based transcript coordinate of the 5'-most nucleotide of
# query_exons (e.g. the start codon of an ORF or CDS), given the transcript's
# exon structure tx_exons.  Splice junctions are handled correctly: only the
# exonic bases are counted, so an ORF and a CDS lying in different exons will
# have their true in-transcript distance reported.
# Returns NA_integer_ if the 5' position falls outside all tx_exons.
.genomic_5prime_to_tx_coord <- function(query_exons, tx_exons, strnd) {
  if (strnd == "+") {
    g         <- min(start(query_exons))
    tx_sorted <- sort(tx_exons)
    in_exon   <- start(tx_sorted) <= g & end(tx_sorted) >= g
    if (!any(in_exon)) return(NA_integer_)
    ex_idx       <- which(in_exon)[1L]
    bases_before <- if (ex_idx > 1L) sum(width(tx_sorted[seq_len(ex_idx - 1L)])) else 0L
    as.integer(bases_before + (g - start(tx_sorted[ex_idx])))
  } else {
    # "-" strand: 5' end is the highest genomic coordinate
    g         <- max(end(query_exons))
    tx_sorted <- rev(sort(tx_exons))   # descending order => first element = 5'-most
    in_exon   <- start(tx_sorted) <= g & end(tx_sorted) >= g
    if (!any(in_exon)) return(NA_integer_)
    ex_idx       <- which(in_exon)[1L]
    bases_before <- if (ex_idx > 1L) sum(width(tx_sorted[seq_len(ex_idx - 1L)])) else 0L
    as.integer(bases_before + (end(tx_sorted[ex_idx]) - g))
  }
}

# Vectorised version of .genomic_5prime_to_tx_coord.
# g_vec  : integer vector of 5' genomic positions
#           (min(start(orf_exons)) for "+", max(end(orf_exons)) for "-")
# tx_exons: GRanges of transcript exons (any order; sorted internally)
# strnd  : single string "+" or "-"
# Returns an integer vector of 0-based transcript coordinates; NA where
# g_vec falls outside all exons.
#
# Uses findInterval() for O(log n_exons) per query instead of a per-element
# R loop, which matters when called for many ORFs on the same transcript.
.genomic_5prime_to_tx_coord_vec <- function(g_vec, tx_exons, strnd) {
  g_vec <- as.integer(g_vec)

  if (strnd == "+") {
    tx_s  <- sort(tx_exons)
    ex_s  <- as.integer(start(tx_s))
    ex_e  <- as.integer(end(tx_s))
    ex_w  <- as.integer(width(tx_s))
    cum_b <- c(0L, cumsum(ex_w)[-length(ex_w)])   # bases before each exon

    # findInterval: k = last exon index where ex_s[k] <= g
    k     <- findInterval(g_vec, ex_s)
    valid <- k >= 1L & k <= length(ex_s) & g_vec <= ex_e[pmax(k, 1L)]
    result <- rep(NA_integer_, length(g_vec))
    result[valid] <- cum_b[k[valid]] + (g_vec[valid] - ex_s[k[valid]])

  } else {
    # "-" strand: tx_sorted descending (5'→3' in transcript coords)
    tx_s  <- rev(sort(tx_exons))
    ex_s  <- as.integer(start(tx_s))
    ex_e  <- as.integer(end(tx_s))    # strictly DECREASING
    ex_w  <- as.integer(width(tx_s))
    cum_b <- c(0L, cumsum(ex_w)[-length(ex_w)])

    # -ex_e is ASCENDING (ex_e is decreasing), so findInterval works directly.
    # k = number of exons where ex_e >= g  (i.e., the exon CONTAINING g is at
    # index k if ex_s[k] <= g).
    k     <- findInterval(-g_vec, -ex_e)
    valid <- k >= 1L & k <= length(ex_s) & ex_s[pmax(k, 1L)] <= g_vec
    result <- rep(NA_integer_, length(g_vec))
    result[valid] <- cum_b[k[valid]] + (ex_e[k[valid]] - g_vec[valid])
  }

  result
}

get_all_orfs <- function(gene_id, annotated_orfs_tab){
  annotated_orfs_tab[annotated_orfs_tab$gene_id == gene_id,]
}

get_orf_isoform_matrices <- function(dtu_results_filtered, annotated_orfs_tab){

  # DTU vector
  dtu <- dtu_results_filtered %>%
    dplyr::select()

}


#### Draft functions for Kozak sequence score calculation
simplify_grl <- function(grl){

  names <- names(grl)

  grl <- map2(grl, name, function(x, y){
    x <- GRanges(seqnames = seqnames(x),
                 IRanges(ranges(x)),
                 strand = strand(x))
    names(x) <- rep(y, length(x))
    return(x)
      }
    )

  grl <- GRangesList(grl)

  return(grl)
}

get_start_position <- function(grl){

  starts <- sapply(grl, function(x){
    min(start(x))
  })

  return(starts)

}

get_stop_position <- function(grl){

  ends <- sapply(grl, function(x){
    max(end(x))
  })

  return(ends)

}
