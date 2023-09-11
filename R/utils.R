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
