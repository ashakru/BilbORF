#' Annotate ORFs with Translation Status Across Transcript Isoforms
#'
#' Maps Open Reading Frames (ORFs) to overlapping transcript isoforms and
#' determines their translation competency by analyzing start codons, stop codons,
#' and internal stops. Each ORF-transcript pair is evaluated independently.
#'
#' @param annotations List containing genomic annotations, as returned by
#'   \code{\link{prepare_annotations_fromGTF}}. Must contain a \code{transcripts}
#'   element with transcript models as a GRangesList.
#' @param orfs GRanges or GRangesList with genomic coordinates of ORF start and
#'   stop codons. For GRanges input, ranges will be grouped by name. Each ORF
#'   can span multiple exons (for spliced ORFs).
#' @param BSgenome A BSgenome object matching the genome build used for annotations.
#'   Required for extracting ORF nucleotide sequences.
#' @param transcript_meta A data.frame with transcript metadata. Must contain
#'   a \code{transcript_id} column matching \code{names(annotations$transcripts)}.
#'   Commonly includes: gene_name, gene_id, transcript_type. This information
#'   is joined to the output table.
#' @param orfs_meta Optional data.frame with ORF metadata. Must contain an
#'   \code{ORF_id} column matching ORF names. Additional columns are joined
#'   to the output.
#' @param start_codons Character vector of valid start codons. Default is
#'   \code{c("ATG","TTG","CTG","GTG")} covering canonical and alternative starts.
#' @param stop_codons Character vector of valid stop codons. Default is
#'   \code{c("TAG", "TAA", "TGA")} covering all standard stop codons.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{table}{A data.frame (tibble) with one row per ORF-transcript pair.
#'       Key columns include:
#'       \itemize{
#'         \item \code{ORF_isoform_id}: Unique identifier for ORF-transcript pair
#'         \item \code{ORF_id}: Original ORF identifier
#'         \item \code{transcript_id}: Overlapping transcript ID
#'         \item \code{gene_id}, \code{gene_name}: Gene annotations
#'         \item \code{seq_nt}: Nucleotide sequence of the ORF in this isoform
#'         \item \code{seq_aa}: Translated amino acid sequence
#'         \item \code{start_codon}, \code{stop_codon}: First and last codons
#'         \item \code{orf_status}: Translation status (see Details)
#'         \item \code{unique_tx_iso}: Groups isoforms producing identical proteins
#'       }}
#'     \item{ranges}{A GRangesList with genomic coordinates of ORF-transcript
#'       intersections, representing the actual translated regions.}
#'   }
#'
#' @details
#' \strong{ORF Status Categories:}
#' \describe{
#'   \item{translatable}{Valid start and stop codons present, no internal stops.
#'     The ORF can produce a full-length protein in this isoform.}
#'   \item{internal_stop}{Contains one or more stop codons before the final
#'     position. Translation would terminate prematurely.}
#'   \item{no_stop}{Missing a valid stop codon at the 3' end. ORF extends to
#'     transcript terminus without proper termination.}
#'   \item{no_start}{Missing a valid start codon at the 5' end but has a stop.
#'     Cannot initiate translation properly.}
#'   \item{no_stop_no_start}{Missing both valid start and stop codons.}
#' }
#'
#' \strong{Algorithm Overview:}
#' \enumerate{
#'   \item Find overlaps between ORFs and transcripts
#'   \item Calculate spliced intersections (ORF coordinates within each transcript)
#'   \item Extract and translate sequences
#'   \item Classify translation status
#'   \item Group isoforms producing identical proteins
#' }
#'
#' @note ORF coordinates should represent the complete ORF boundaries (start to
#'   stop codon), not just the coding sequence. For multi-exon ORFs, provide
#'   all exonic segments.
#'
#' @import ORFik
#' @import tibble
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom dplyr mutate case_when left_join group_by relocate
#' @export
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(rtracklayer)
#' 
#' # Prepare annotations
#' gtf <- "path/to/annotation.gtf"
#' annotations <- prepare_annotations_fromGTF(gtf, BSgenome.Hsapiens.UCSC.hg38)
#' 
#' # Load ORFs from BED file
#' orfs <- import("path/to/orfs.bed", format = "BED")
#' names(orfs) <- orfs$name
#' 
#' # Prepare transcript metadata
#' transcripts_meta <- import(gtf, format = "GTF") %>%
#'   as.data.frame() %>%
#'   filter(type == "transcript") %>%
#'   select(gene_name, gene_id, transcript_id, transcript_type) %>%
#'   distinct()
#' 
#' # Annotate ORFs
#' result <- annotate_orf_isoforms(
#'   annotations, 
#'   orfs, 
#'   BSgenome.Hsapiens.UCSC.hg38,
#'   transcripts_meta
#' )
#' 
#' # Explore results
#' table(result$table$orf_status)
#' head(result$table)
#' }
#' 
#' # Using package example data
#' gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", 
#'                    package = "BilbORF")
#' bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
#' 
#' if (gtf != "" && bed != "") {
#'   library(BSgenome.Hsapiens.UCSC.hg38)
#'   library(rtracklayer)
#'   
#'   BSgenome <- BSgenome.Hsapiens.UCSC.hg38
#'   annotations <- prepare_annotations_fromGTF(gtf, BSgenome)
#'   orfs <- import(bed, format = "BED")
#'   names(orfs) <- orfs$name
#'   
#'   transcripts_meta <- import(gtf, format = "GTF") %>%
#'     as.data.frame() %>%
#'     dplyr::filter(type == "transcript") %>%
#'     dplyr::select(gene_name, gene_id, transcript_id, transcript_type) %>%
#'     dplyr::distinct()
#'   
#'   annotated_orfs <- annotate_orf_isoforms(
#'     annotations, orfs[1:5], BSgenome, transcripts_meta
#'   )
#'   
#'   print(table(annotated_orfs$table$orf_status))
#' }
#'
#' @seealso 
#' \code{\link{prepare_annotations_fromGTF}} for creating annotation input
#' \code{\link{diff_orf_usage}} for differential usage analysis
annotate_orf_isoforms <- function(annotations, orfs, BSgenome,
                                  transcript_meta, orfs_meta = NULL,
                                  start_codons = c("ATG","TTG","CTG","GTG"),
                                  stop_codons = c("TAG", "TAA", "TGA")){

  tx <- annotations$transcripts

  # Group ORFs by names
  if (class(orfs) != "GRangesList"){
    orfs <- ORFik::groupGRangesBy(orfs)
  }

  # Link ORFs with overlapping transcripts
  orf_tx_overlap <- findOverlaps(orfs, tx)

  # Get spliced coordinates
  orf_by_tx <- orfs[queryHits(orf_tx_overlap)]
  tx_by_orf <- tx[subjectHits(orf_tx_overlap)]
  orf_in_tx <- GenomicRanges::intersect(tx_by_orf, orf_by_tx)

  # Rank exons
  orf_in_tx <- rank_exons(orf_in_tx)

  # Set ORF names
  ORF_isoform_id <- names(orf_in_tx) <- paste(names(orf_by_tx), paste0("iso-", names(tx_by_orf)), sep = "_")
  names(orf_in_tx) <- ORF_isoform_id

  # Get AA seqeunces
  seq_nt <- extractTranscriptSeqs(BSgenome, orf_in_tx)
  seq_aa <- suppressWarnings(translate(seq_nt))

  # Computing Kozak sequence (currently not working)
  # orf_test <- simplify_grl(orf_in_tx[1:2])
  # orf_test <- makeORFNames(orf_in_tx[1:2])
  # names(orf_test) <- names(tx_by_orf[1:2])
  # tx_test <- tx_by_orf[names(orf_test)]
  # tx_test <- makeORFNames(tx_test)
  #
  # kozakSequenceScore(orf_test[1], tx_test[[1]], BSgenome)
  # pmapToTranscriptF(orf_test, tx[names(orf_test)])
  # sequences <- startRegionString(orf_test, tx[names(orf_test)], BSgenome, 2, 3)
  # grl <- startRegion(orf_test, tx[names(orf_test)], is.sorted = TRUE, 9, 5)

  # Assign ORF status
  orf_status <- tibble(ORF_isoform_id = names(orf_in_tx),
                       ORF_id = names(orf_by_tx),
                       transcript_id = names(tx_by_orf),
                       seq_nt = as.character(seq_nt),
                       seq_aa = as.character(seq_aa),
                       len_nt = nchar(seq_nt),
                       len_aa = nchar(gsub("*", "", seq_aa)),
                       seqnames = sapply(orf_in_tx, function(x){unique(seqnames(x))}),
                       start = get_start_position(orf_in_tx),
                       end = get_stop_position(orf_in_tx)) %>%
    mutate(length_nt = nchar(seq_nt),
           length_aa = nchar(gsub("*", "", seq_aa)),
           start_codon = substr(seq_nt, 1, 3),
           stop_codon = substr(seq_nt, length_nt-2, length_nt),
           last_frame = length_nt %% 3,
           orf_status = case_when(grepl("[*]", substr(seq_aa, 1, length_aa-1)) ~ "internal_stop",
                                  !(start_codon %in% start_codons) & !(stop_codon %in% stop_codons) ~ "no_stop_no_start",
                                  !(stop_codon %in% stop_codons) ~ "no_stop",
                                  !(start_codon %in% start_codons) ~ "no_start",
                                  (start_codon %in% start_codons) & (stop_codon %in% stop_codons) ~ "translatable")) %>%
    group_by(ORF_id, orf_status) %>%
    mutate(unique_iso_n =  as.numeric(factor(seq_aa)),
           unique_tx_iso = paste(orf_status, unique_iso_n, sep = "_"))

  # Join transcript metadata table, if present
 orf_status <- orf_status %>%
      left_join(transcript_meta, by = "transcript_id")


  # Join ORF metadata table, if present
  if (!is.null(orfs_meta)) {
    orf_status <- orf_status %>%
      left_join(orf_meta, by = "ORF_id")
  }

 # Relocate the most important columns
 orf_status <- orf_status %>%
   dplyr::relocate(ORF_isoform_id, ORF_id, transcript_id, gene_id, gene_name,
                   seqnames, start, end, start_codon, stop_codon, orf_status)

  return(list(table = orf_status,
              ranges = orf_in_tx))

}


