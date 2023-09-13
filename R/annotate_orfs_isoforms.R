#' Annotate ORF with their status respective of
#'
#' @param annotations List with genomic annotations, output of `prepare_annotations()`
#' @param orfs GRanges or GRangesList with genomic coordinates of ORFs of interest
#' @param transcript_meta (optional) A data.frame with transcript metadata (for example: gene_name, gene_id, transcript_type), must contain a column `transcript_id` compatible with `names(annotations$transcripts)`
#' @param orfs_meta (optional) A data.frame with ORF metadata, must contain a column `ORF_id` compatible with `names(orfs)`
#' @param start_codons A vector with codons recognized as start codon; default `c("ATG","TTG","CTG","GTG")`
#' @param stop_codons A vector with codons recognized as stop codon; default `c("TAG", "TAA", "TGA")`
#'
#' @return A list containing a table with annotated ORF-transcript pairs (`$table`) and GRangesList with genomic coordinates of identified ORFs (`$ranges`)
#' @import GenomicRanges
#' @import ORFik
#' @export
#'
#' @examples
#' BSgenome <- BSgenome.Hsapiens.UCSC.hg38
#' gtf <- "inst/extdata/gencode.v35.annotation_chr10.gtf"
#' bed <- "inst/extdata/Ribo-seq_ORFs.bed"
#'
#' # Prepare annotations
#' annotations <- prepare_annotations_fromGTF(gtf, BSgenome)
#'
#' # Pull transcripts metadata
#' transcripts_meta <- import(gtf, format = "GTF") %>%
#'   as.data.frame() %>%
#'   dplyr::filter(type == "transcript") %>%
#'   dplyr::select(gene_name, transcript_id, transcript_type) %>%
#'   distinct()
#'
#' # Load ORFs
#' orfs <- import(bed, format = "BED")
#' names(orfs) <- orfs$name
#' orf_tab <- as.data.frame(orfs)
#' annotated_orfs <- annotate_orf_isoforms(annotations, orfs, BSgenome, transcripts_meta)
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

  # Computing Kozak sequence
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
                                  (start_codon %in% start_codons) & (stop_codon %in% stop_codons) ~ "translatable"))

  # Join transcript metadata table, if present
 orf_status <- orf_status %>%
      left_join(transcript_meta, by = "transcript_id")


  # Join ORF metadata table, if present
  if (!is.null(orfs_meta)) {
    orf_status <- orf_status %>%
      left_join(orf_meta, by = "ORF_id")
  }

  return(list(table = orf_status,
              ranges = orf_in_tx))

}


