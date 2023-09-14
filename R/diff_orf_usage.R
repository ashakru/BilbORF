#' Calculate Differential ORF Usage
#'
#' This function calculates differential ORF usage based on the provided inputs.
#'
#' @param annotated_orfs_tab A data frame containing information about annotated ORFs.
#' @param dtu_results A data frame containing results of differential transcript usage analysis.
#' @param selected_genes A character vector of gene IDs for which differential ORF usage will be calculated.
#' @param orf_translation (optional) A data frame containing ORF translations.
#'
#' @return A list containing two elements: 'gene_stats' and 'gene_tables'.
#'   - 'gene_stats': A data frame with statistics about the differential ORF usage for each gene.
#'   - 'gene_tables': A list of data tables, each containing the ORF/isoform combinations for a specific gene.
#'
#' @examples
#' # Example Usage
#' diff_orf_usage(annotated_orfs_tab, dtu_results, selected_genes)
#'
#' @export
#'
#' @import data.table
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
