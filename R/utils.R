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

