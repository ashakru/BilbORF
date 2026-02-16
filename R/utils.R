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

    grl <- lapply(grl, function(x){
      str <- as.character(strand(x))
      if(all(str == "-")) {
        x$exon_rank <- length(x):1
      } else {
        x$exon_rank <- 1:length(x)
      }

      x <- x[order(x$exon_rank, decreasing = F)]

      return(x)
    })

  grl <- GRangesList(grl)

  return(grl)

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
