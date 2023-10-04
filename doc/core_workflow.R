## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----out.width = '70%', echo = FALSE------------------------------------------
knitr::include_graphics("img/orf_translation_splicing.png")

## ----setup, message=F---------------------------------------------------------
library(BilbORF)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)

## ----annotations--------------------------------------------------------------
# Load reference
BSgenome <- BSgenome.Hsapiens.UCSC.hg38
gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")

# Prepare annotations
annotations <- prepare_annotations_fromGTF(gtf, BSgenome)

## ----orf_status---------------------------------------------------------------
# Load ORFs
bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
orfs <- import(bed, format = "BED")
names(orfs) <- orfs$name
orf_tab <- as.data.frame(orfs)

# Pull transcripts metadata
transcripts_meta <- import(gtf, format = "GTF") %>%
  as.data.frame() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(gene_name, gene_id, transcript_id, transcript_type) %>%
  dplyr::distinct()


# Annotate ORF status on all overlapping transcripts 
annotated_orfs <- annotate_orf_isoforms(annotations, orfs, BSgenome, transcripts_meta)

## ----message=FALSE------------------------------------------------------------
# Get an example BAM file with Ribo-Seq reads
bam <- system.file("extdata", "riboseq_chr10.bam", package = "BilbORF")

# Define offsets
offsets <- data.frame(fraction = 25:30,
                      offsets_start = -12) 

# Quantify ORFs
orf_translation <- count_p_sites(bam, offsets, annotated_orfs)

## ----message=FALSE------------------------------------------------------------
library(DRIMSeq)

# Pull transcript metadata
transcripts_table <- read.csv(system.file("extdata", "gencode.v35.tx2gene.csv", package = "BilbORF"))
transcripts_table <- transcripts_table %>%
  dplyr::mutate(transcript_ensembl_id = sapply(strsplit(transcript_id, "[.]"), function(x){x[1]}))

# Differential transcript usage
d <- dmDSdata(counts = gtex_data$counts, samples = gtex_data$meta)
d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 10)
design_full <- model.matrix(~ group, data = samples(d))

set.seed(123)
d <- dmPrecision(d, design = design_full)
d <- dmFit(d, design = design_full, verbose = 1)
d <- dmTest(d, coef = "groupLung", verbose = 1)

# Tidy results
dtu_results <- results(d, level = "feature") %>%
  left_join(proportions(d)) %>%
  dplyr::mutate(mean_lung = rowMeans(dplyr::select(., starts_with("Lung")), na.rm = TRUE),
         mean_liver = rowMeans(dplyr::select(., starts_with("Liver")), na.rm = TRUE),
         log2FoldChange = log2(mean_lung/mean_liver)) %>%
  dplyr::rename(transcript_id = feature_id) %>%
  dplyr::arrange(adj_pvalue) %>%
  dplyr::left_join(dplyr::select(transcripts_table, transcript_ensembl_id, gene_name), 
            by = c("transcript_id" = "transcript_ensembl_id")) %>%
  dplyr::relocate(gene_name, gene_id, transcript_id, pvalue, adj_pvalue, mean_lung, mean_liver, log2FoldChange) 

## ----message=FALSE------------------------------------------------------------
dtu_results <- dtu_results %>%
  mutate(rank = log2FoldChange)

annotated_orfs_tab <- annotated_orfs$table %>%
  mutate(transcript_id = sapply(strsplit(transcript_id, "[.]"), 
                                            function(x){x[1]}))
  
selected_genes <- dtu_results %>%
  mutate(hasORF = transcript_id %in% annotated_orfs_tab$transcript_id) %>%
  dplyr::filter(hasORF & adj_pvalue < 0.05)
selected_genes <- unique(selected_genes$gene_id)

orf_landscape <- diff_orf_usage(annotated_orfs_tab, dtu_results, selected_genes)

## ----message=FALSE------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# Selected genes heatmap 
gene_stats <- orf_landscape$gene_stats
mat_unscaled <- as.tibble(gene_stats) %>%
  dplyr::select(-gene_id) %>%
  column_to_rownames("gene_name")
mat <- scale(mat_unscaled)

# Graphics
col_fun <- colorRamp2(c(-2, 0, 2), rev(c("#4D004B", "#8C96C6", "#9EBCDA")))
col_fun(seq(-3, 3))

# Heatmap
Heatmap(mat,
        cluster_rows = F,
        cluster_columns = F,
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(round(mat_unscaled[i, j], 2), x, y, gp = gpar(col = "white")) # Set text color to white
        })

# Examples: NUTM2B-AS1 and ADK
examples <- c("ENSG00000225484.7", "ENSG00000156110.14")
col_dtu <-  colorRamp2(c(-2, 0, 2), c("#225EA8", "white", "#E31A1C"))
col_orf <-  structure(c("#4D004B", "#8C6BB1", "#D4B9DA", "#8C96C6", "#9EBCDA","#BFD3E6", "grey"),
                    names = c("translatable_1", "translatable_2", "internal_stop", "no_stop_no_start", "no_stop", "orf_loss"))

for (e in seq_along(examples)){
  gene_tab <- orf_landscape$gene_tables[[examples[e]]]
  gene_name <- gene_stats$gene_name[gene_stats$gene_id == examples[e]]
  
  dtu_results <- as.tibble(gene_tab) %>%
    dplyr::select(transcript_id, rank) %>%
    column_to_rownames("transcript_id")
  dtu_results <- as.matrix(dtu_results)
  
  orf_tab <- as.tibble(gene_tab) %>%
    column_to_rownames("transcript_id") %>%
    dplyr::select(-adj_pvalue, -rank) 
  orf_tab <- as.matrix(orf_tab)  
  
  dtu_heatmap <- Heatmap(dtu_results, 
               name = "DTU",
               col = col_dtu,
               cluster_columns = F,
               cluster_rows = F,
               row_title_gp = gpar(fontsize = 8),
               show_row_names = T, 
               row_names_side = "left",  
               border = TRUE,
               border_gp = gpar(col = "#838383", lwd = 0.5))
  
  orf <- Heatmap(orf_tab, 
               name = "ORFs",
               col = col_orf, 
               cluster_rows = F, 
               cluster_columns = T,
               row_title_gp = gpar(fontsize = 8),
               show_row_names = F,
               rect_gp = gpar(col = "white", lwd = 1),
               height = unit(3, "cm"))
  
  ht_list <- dtu_heatmap + orf
  
  draw(ht_list, column_title = gene_name)
  
}



## -----------------------------------------------------------------------------
sessionInfo()

