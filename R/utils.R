#' Title
#'
#' @param x,y The object from/on which to get/set the seqlevels and seqlevels style
#'
#' @return NULL
#'
#' @examples check_seq_levels(Granges_1, Granges_2)
#' @keywords internal
#' @export
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

#' Function to assign appropriate exon ranking depending on the ranges strand
#'
#' @param grl GRangesList object
#' @return A list of genomic ranges objects with exon ranks assigned.
#' @export
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
