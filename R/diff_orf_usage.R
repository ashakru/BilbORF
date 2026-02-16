#' Analyze Differential ORF Usage from Transcript Isoform Changes
#'
#' Links differential transcript usage (DTU) results to changes in ORF translation
#' status. Creates matrices showing which ORFs are gained, lost, or switched
#' between translatable and non-translatable states across isoform changes.
#'
#' @param annotated_orfs_tab Data frame from \code{annotate_orf_isoforms()$table}
#'   containing ORF-transcript annotations. Must include columns: ORF_id,
#'   transcript_id, gene_id, unique_tx_iso.
#' @param dtu_results Data frame with differential transcript usage results.
#'   Required columns:
#'   \describe{
#'     \item{transcript_id}{Transcript identifiers matching annotated_orfs_tab}
#'     \item{gene_id}{Gene identifiers}
#'     \item{gene_name}{Gene symbols for labeling}
#'     \item{pvalue}{Statistical significance of DTU}
#'     \item{adj_pvalue}{Adjusted p-value (e.g., from multiple testing correction)}
#'     \item{rank}{Numeric ranking metric (e.g., log2FoldChange or ΔUsage)}
#'   }
#'   Typically from DRIMSeq, DEXSeq, or similar DTU tools.
#' @param selected_genes Character vector of gene IDs to analyze. Usually
#'   genes with significant DTU events and detectable ORFs.
#' @param orf_translation Optional data frame with ORF translation metrics
#'   from \code{\link{count_p_sites}}. Currently not used in calculations
#'   but reserved for future quantitative analysis.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{gene_stats}{Data frame with one row per gene, sorted by n_events
#'       (descending). Columns:
#'       \itemize{
#'         \item gene_id, gene_name: Gene identifiers
#'         \item n_tx_isoforms: Number of transcript isoforms for this gene
#'         \item n_signif_dtu: Count of isoforms with adj_pvalue < 0.05
#'         \item n_orfs: Number of distinct ORFs overlapping this gene
#'         \item n_events: Total ORF-isoform combinations (potential translation events)
#'         \item n_translatable: Count of translatable ORF-isoform pairs
#'       }}
#'     \item{gene_tables}{Named list of data.tables, one per gene. Each table
#'       is a matrix with:
#'       \itemize{
#'         \item Rows: Transcript isoforms (with adj_pvalue, rank metadata)
#'         \item Columns: ORFs
#'         \item Values: Translation status (e.g., "translatable_1", "internal_stop",
#'           "orf_loss" if ORF doesn't overlap that isoform)
#'       }}
#'   }
#'
#' @details
#' \strong{Analysis Strategy:}
#' For each gene, the function constructs a matrix showing the fate of each ORF
#' across all transcript isoforms. This reveals:
#' \itemize{
#'   \item \strong{ORF gain/loss}: ORFs present in some isoforms but not others
#'   \item \strong{Status switches}: ORFs changing between translatable and
#'     non-translatable due to isoform-specific sequences
#'   \item \strong{Protein isoform diversity}: Multiple translatable ORF states
#'     producing different protein products
#' }
#'
#' The \code{unique_tx_iso} column groups ORF-isoform pairs producing identical
#' amino acid sequences, helping identify truly distinct protein products.
#'
#' \strong{Typical Workflow:}
#' \enumerate{
#'   \item Perform DTU analysis (e.g., with DRIMSeq)
#'   \item Annotate ORFs with \code{annotate_orf_isoforms()}
#'   \item Filter to genes with significant DTU and detected ORFs
#'   \item Run \code{diff_orf_usage()} to link isoform changes to ORF changes
#'   \item Visualize results with heatmaps (see vignette)
#' }
#'
#' @note This function provides qualitative analysis of ORF status changes.
#'   Quantitative changes in ORF abundance require combining with Ribo-seq
#'   data and isoform quantification.
#'
#' @export
#' @import data.table
#'
#' @examples
#' \dontrun{
#' # Typical workflow
#' library(DRIMSeq)
#' library(BilbORF)
#' 
#' # 1. Perform DTU analysis
#' d <- dmDSdata(counts = counts, samples = samples)
#' d <- dmFilter(d)
#' d <- dmPrecision(d, design = design)
#' d <- dmFit(d, design = design)
#' d <- dmTest(d, coef = "condition")
#' 
#' # 2. Format DTU results
#' dtu_results <- results(d, level = "feature") %>%
#'   mutate(rank = log2FoldChange) # or another ranking metric
#' 
#' # 3. Annotate ORFs (see ?annotate_orf_isoforms)
#' annotated_orfs <- annotate_orf_isoforms(...)
#' 
#' # 4. Select genes with significant DTU and ORFs
#' selected_genes <- dtu_results %>%
#'   filter(adj_pvalue < 0.05,
#'          transcript_id %in% annotated_orfs$table$transcript_id) %>%
#'   pull(gene_id) %>%
#'   unique()
#' 
#' # 5. Analyze differential ORF usage
#' orf_landscape <- diff_orf_usage(
#'   annotated_orfs$table,
#'   dtu_results,
#'   selected_genes
#' )
#' 
#' # 6. Explore results
#' head(orf_landscape$gene_stats)
#' orf_landscape$gene_tables[[1]]  # Matrix for first gene
#' }
#' 
#' # Simple example with mock data
#' orf_tab <- data.frame(
#'   ORF_id = c("orf1", "orf1", "orf2"),
#'   transcript_id = c("tx1", "tx2", "tx1"),
#'   gene_id = c("gene1", "gene1", "gene1"),
#'   unique_tx_iso = c("translatable_1", "internal_stop", "translatable_1")
#' )
#' 
#' dtu_results <- data.frame(
#'   transcript_id = c("tx1", "tx2"),
#'   gene_id = c("gene1", "gene1"),
#'   gene_name = c("GENEA", "GENEA"),
#'   pvalue = c(0.01, 0.05),
#'   adj_pvalue = c(0.02, 0.08),
#'   rank = c(-1.5, 2.0)
#' )
#' 
#' results <- diff_orf_usage(orf_tab, dtu_results, "gene1")
#' print(results$gene_stats)
#' print(results$gene_tables$gene1)
#' 
#' @seealso
#' \code{\link{annotate_orf_isoforms}} for preparing ORF input
#' \code{DRIMSeq}, \code{DEXSeq} for DTU analysis
#' Vignette \code{vignette("core_workflow", package = "BilbORF")} for visualization
diff_orf_usage <- function(annotated_orfs_tab, dtu_results, selected_genes, orf_translation = NULL) {

    # Check if dtu_results tab contains all information required
    required_columns <-  c("transcript_id", "gene_id", "gene_name", "pvalue", "rank", "adj_pvalue")

    if (!all(required_columns %in% colnames(dtu_results))) {
      missing <- colnames(dtu_results)[!colnames(dtu_results) %in% required_columns]
      stop(
        paste(
          "Not all required columns are present in dtu_results. Missing columns: ",
          missing
        )
      )
    }

    # Subset set of genes
    dtu_results_filtered <- dtu_results %>%
      dplyr::select(all_of(required_columns)) %>%
      dplyr::filter(gene_id %in% selected_genes) %>%
      dplyr::arrange(gene_id)

    # Gene to name
    gene2name <- dtu_results_filtered %>%
      dplyr::select(gene_id, gene_name) %>%
      dplyr::distinct()

    # Build matrix of ORF-transcript status
    orf_tx_tab <- annotated_orfs_tab %>%
      dplyr::filter(gene_id %in% selected_genes) %>%
      dplyr::arrange(gene_id, transcript_id, ORF_id) %>%
      dplyr::select(ORF_id, transcript_id, unique_tx_iso)
    orf_tx_tab <- data.table::as.data.table(orf_tx_tab)

    orf_tx_tables <- list()
    orf_stats_per_gene <- data.table()
    # Build per gene matrices of ORF/isoform combinations
    for (g in selected_genes){
      gene_name <- gene2name$gene_name[gene2name$gene_id == g]
      iso_tab <- dtu_results_filtered[dtu_results_filtered$gene_id == g,]
      g_iso <- iso_tab$transcript_id
      g_orfs <- na.omit(unique(orf_tx_tab$ORF_id[match(g_iso, orf_tx_tab$transcript_id)]))
      g_mat <- as.data.table(matrix(NA, nrow = length(g_iso), ncol = length(g_orfs)))
      colnames(g_mat) <- g_orfs
      g_mat <- cbind(iso_tab[,c("transcript_id", "adj_pvalue", "rank")], g_mat)

      for (o in g_orfs){
        g_mat[,o] <- orf_tx_tab[match(g_iso, orf_tx_tab$transcript_id), "unique_tx_iso"]
      }

      g_mat[is.na(g_mat)] <- "orf_loss"

      # Return ORF landscape stats per gene
      unique_events <- as.matrix(unique(g_mat[,-1:-3]))

      n_translatable <- sum(grepl("translatable", unique_events))
      n_events <- ncol(unique_events)*nrow(unique_events)

      n_signif_dtu <- sum(iso_tab$adj_pvalue < 0.05)

      gene_tab <- data.table(gene_id = g,
                             gene_name = gene_name,
                             n_tx_isoforms = length(g_iso),
                             n_signif_dtu = n_signif_dtu,
                             n_orfs = length(g_orfs),
                             n_events = n_events,
                             n_translatable = n_translatable)

      orf_stats_per_gene <- rbind(orf_stats_per_gene, gene_tab)

      # Save full matrices
      orf_tx_tables[[g]] <- g_mat
    }

    orf_stats_per_gene <- orf_stats_per_gene[order(orf_stats_per_gene$n_events, decreasing = T),]
    orf_tx_tables <- orf_tx_tables[orf_stats_per_gene$gene_id]

    return(list(gene_stats = orf_stats_per_gene,
                gene_tables = orf_tx_tables))
  }
