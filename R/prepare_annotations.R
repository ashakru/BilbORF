#' Prepare Transcript Annotations from GTF File
#'
#' This function extracts transcript isoform models from a GTF file and returns
#' them as genomic coordinates grouped by transcript. It creates a TxDb object
#' internally and validates sequence level compatibility with the provided genome.
#'
#' @param gtf Character string specifying the path to a GTF file containing
#'   transcript isoform models. The GTF should follow standard format with
#'   exon features annotated by transcript_id.
#' @param BSgenome A BSgenome object compatible with the GTF file. The
#'   seqlevels style (e.g., "chr1" vs "1") and sequence names must match
#'   between the GTF and BSgenome.
#'
#' @return A list with one element:
#'   \item{transcripts}{A GRangesList where each element represents a transcript
#'   and contains all its exons as genomic ranges. Names correspond to transcript IDs.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Builds a TxDb object from the GTF using \code{makeTxDbFromGFF}
#'   \item Validates that sequence levels match between TxDb and BSgenome
#'   \item Extracts exons grouped by transcript using \code{exonsBy}
#' }
#'
#' Sequence level validation ensures chromosome naming consistency. For example,
#' UCSC-style genomes use "chr1" while Ensembl uses "1". Both GTF and BSgenome
#' must use the same convention.
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @export
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' 
#' # Using a local GTF file
#' gtf <- "path/to/gencode.v35.annotation.gtf"
#' annotations <- prepare_annotations_fromGTF(gtf, BSgenome.Hsapiens.UCSC.hg38)
#' 
#' # Access transcript models
#' head(names(annotations$transcripts))  # Transcript IDs
#' annotations$transcripts[[1]]  # First transcript's exons
#' }
#' 
#' # Using package example data
#' gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", 
#'                    package = "BilbORF")
#' if (gtf != "") {
#'   library(BSgenome.Hsapiens.UCSC.hg38)
#'   annotations <- prepare_annotations_fromGTF(gtf, BSgenome.Hsapiens.UCSC.hg38)
#'   print(length(annotations$transcripts))
#' }
#'
#' @seealso \code{\link{annotate_orf_isoforms}} for using these annotations
#'   to map ORFs to transcripts
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
