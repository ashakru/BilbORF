#' Quantify Ribosome P-sites Within ORFs from Ribo-seq Data
#'
#' Analyzes ribosome profiling (Ribo-seq) data to quantify P-site coverage within
#' annotated ORFs. Calculates frame-specific coverage, entropy, and ORF scores
#' to assess translation activity and frame preference.
#'
#' @param bam Character string specifying the path to a BAM file containing
#'   aligned ribosome-protected fragments (RPFs) from Ribo-seq. Must be indexed
#'   (BAM.bai file present).
#' @param offsets A data.frame with P-site offset information. Must contain:
#'   \describe{
#'     \item{fraction}{Read lengths to include (e.g., 25:30 for 25-30nt reads)}
#'     \item{offsets_start}{Numeric offset from read 5' end to P-site position
#'       (typically negative, e.g., -12). Can be estimated using
#'       \code{ORFik::detectRibosomeShifts}}
#'   }
#' @param annotated_orfs List returned by \code{\link{annotate_orf_isoforms}},
#'   containing \code{$table} and \code{$ranges} elements with ORF annotations.
#'
#' @return A data.frame with one row per ORF containing:
#'   \describe{
#'     \item{ORF_id}{ORF identifier}
#'     \item{entropy}{Shannon entropy of read distribution across the three frames.
#'       Values range from 0 (all reads in one frame) to log2(3) ≈ 1.58 (uniform
#'       distribution). Lower entropy indicates stronger frame preference.}
#'     \item{frame_0, frame_1, frame_2}{Number of P-sites in each reading frame.
#'       Frame 0 is the expected translation frame.}
#'     \item{frames_sum}{Total P-sites across all frames}
#'     \item{orfScore}{ORF score measuring translation confidence based on
#'       frame preference (only for ORFs ≥3nt)}
#'   }
#'
#' @details
#' \strong{Analysis Steps:}
#' \enumerate{
#'   \item Load ribosome footprints from BAM
#'   \item Filter reads by specified lengths (offsets$fraction)
#'   \item Calculate P-site positions using provided offsets
#'   \item Collapse duplicate reads (summing their counts)
#'   \item Compute frame-specific coverage for each ORF
#'   \item Calculate entropy and ORF scores
#' }
#'
#' \strong{Interpreting Results:}
#' \itemize{
#'   \item High frame_0 coverage with low entropy suggests active translation
#'   \item Uniform frame distribution (high entropy) may indicate noise or
#'     non-specific binding
#'   \item ORF score combines frame preference with coverage for translation
#'     evidence
#' }
#'
#' @note
#' \itemize{
#'   \item P-site offsets are critical for accurate positioning. Use
#'     \code{ORFik::detectRibosomeShifts} to estimate offsets for your data.
#'   \item Results are aggregated across all ORF isoforms (not isoform-specific).
#'   \item ORFs shorter than 3nt are excluded from orfScore calculation.
#' }
#'
#' @import ORFik
#' @importFrom dplyr rename summarise
#' @export
#'
#' @examples
#' \dontrun{
#' # Typical workflow
#' library(BilbORF)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' 
#' # 1. Prepare annotations and annotate ORFs (see ?annotate_orf_isoforms)
#' annotations <- prepare_annotations_fromGTF(gtf, BSgenome.Hsapiens.UCSC.hg38)
#' annotated_orfs <- annotate_orf_isoforms(annotations, orfs, 
#'                                         BSgenome.Hsapiens.UCSC.hg38,
#'                                         transcripts_meta)
#' 
#' # 2. Define P-site offsets for your data
#' offsets <- data.frame(
#'   fraction = 25:30,      # Read lengths to include
#'   offsets_start = -12    # Offset to P-site
#' )
#' 
#' # 3. Quantify P-sites
#' orf_translation <- count_p_sites(
#'   "path/to/riboseq.bam",
#'   offsets,
#'   annotated_orfs
#' )
#' 
#' # 4. Analyze results
#' head(orf_translation)
#' hist(orf_translation$entropy)
#' plot(orf_translation$frames_sum, orf_translation$entropy)
#' }
#' 
#' # Using package example data
#' bam <- system.file("extdata", "riboseq_chr10.bam", package = "BilbORF")
#' gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", 
#'                    package = "BilbORF")
#' bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
#' 
#' if (bam != "" && gtf != "" && bed != "") {
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
#'     annotations, orfs[1:20], BSgenome, transcripts_meta
#'   )
#'   
#'   offsets <- data.frame(fraction = 25:30, offsets_start = -12)
#'   results <- count_p_sites(bam, offsets, annotated_orfs)
#'   print(head(results))
#' }
#' 
#' @seealso 
#' \code{\link{annotate_orf_isoforms}} for preparing ORF annotations
#' \code{ORFik::detectRibosomeShifts} for estimating P-site offsets
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
    dplyr::summarise(p_sites = sum(score)) %>%
    tidyr::spread(frame, p_sites) %>%
    dplyr::rename(ORF_id = genes) %>%
    mutate(frames_sum = frame_0 + frame_1 + frame_2)

  # Concentrate results table
  results_table <- entropy_per_orf %>%
    left_join(cov_per_frame) %>%
    left_join(orf_score)

  return(results_table)
}
