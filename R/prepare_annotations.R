#' Prepare annotation set
#'
#' @param gtf Path to a GTF file including transcript isoform models
#' @param BSgenome BSgenome object compatible with provided GTF file
#'
#' @return GRangesList with genomic coordinates of different genomic regions
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @export
#'
#' @examples
#'
#' library(BSgenome.Hsapiens.UCSC.hg38) # Load selected BSgenome
#' gtf <- "gencode.v29.annotations.gtf"
#' annotations <- prepare_annotations_fromGTF(gtf, BSgenome.Hsapiens.UCSC.hg38)
prepare_annotations_fromGTF <- function(gtf, BSgenome) {

  # Build transcripts models
  txdb <- suppressWarnings(suppressMessages(makeTxDbFromGFF(gtf)))

  check_seq_levels(txdb, BSgenome)
  transcripts <- exonsBy(txdb, by = "tx", use.names=T)
  annotations <- list(transcripts = transcripts)

  return(annotations)

}


# get_transcript_biotypes <- function(gtf){
#
#   tx2biotypes  <- import(gtf, format = "GTF") %>%
#     as.data.frame() %>%
#     dplyr::filter(type == "transcript") %>%
#     dplyr::select(gene_name, transcript_id, transcript_type) %>%
#     distinct()
#
#   return(tx2biotypes)
#
# }
