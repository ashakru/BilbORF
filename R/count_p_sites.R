#' Count P-sites within selected ORFs
#'
#' This function analyzes sequencing data stored in a BAM file format to identify P-sites (ribosome profiling sites) within annotated open reading frames (ORFs).
#'
#' @param bam Path to the BAM file containing the sequencing data.
#' @param offsets A data frame containing information about the fractional offsets for P-site identification.
#' @param annotated_orfs A data frame containing genomic coordinates of annotated ORFs.
#'
#' @return A data frame containing computed metrics related to P-sites, entropy, and ORF scores per frame.
#' @import ORFik
#' @examples
#' results <- count_p_sites("example.bam", offsets_df, annotated_orfs_df)
#' @export
#'
count_p_sites <- function(bam, offsets, annotated_orfs){

  # Load BAM
  footprints <- readBam(bam)

  # Find P-sites
  footprints <- footprints[readWidths(footprints) %in% as.numeric(offsets$fraction)]
  p_sites <- shiftFootprints(footprints, offsets)
  p_sites <- collapseDuplicatedReads(p_sites, addSizeColumn = TRUE)


  # Compute read coverage entropy per frame
  orf_grl <- annotated_orfs$ranges
  entropy_per_orf <- data.frame(ORF_id = names(orf_grl),
                                entropy = entropy(orf_grl, p_sites, weight = "score"))

  orf_grl_filtered <- orf_grl[sum(width(orf_grl)) >= 3]

  # Compute ORF score per frame
  orf_score <- orfScore(orf_grl_filtered, p_sites)
  orf_score$ORF_id <- names(orf_grl_filtered)
  orf_score <- orf_score[,-1:-3]

  # Compute p-sites coverage per frame
  cov_per_frame <- regionPerReadLength(orf_grl,
                                     p_sites,
                                     withFrames = T,
                                     scoring = "frameSumPerLG",
                                     exclude.zero.cov.grl = T) %>%
    group_by(genes, frame) %>%
    mutate(frame = paste0("frame_", frame)) %>%
    summarise(p_sites = sum(score)) %>%
    tidyr::spread(frame, p_sites) %>%
    rename(ORF_id = genes) %>%
    mutate(frames_sum = frame_0 + frame_1 + frame_2)

  # Concentrate results table
  results_table <- entropy_per_orf %>%
    left_join(cov_per_frame) %>%
    left_join(orf_score)

  return(results_table)
}
