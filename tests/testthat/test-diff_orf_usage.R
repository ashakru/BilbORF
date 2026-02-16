test_that("diff_orf_usage calculates differential usage correctly", {
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
  skip_if(gtf == "" || bed == "", "Test data not found")
  
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(rtracklayer))
  BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  
  # Prepare test data
  annotations <- suppressWarnings(suppressMessages(
    prepare_annotations_fromGTF(gtf, BSgenome)
  ))
  
  orfs <- import(bed, format = "BED")
  names(orfs) <- orfs$name
  
  transcripts_meta <- import(gtf, format = "GTF") %>%
    as.data.frame() %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::select(gene_name, gene_id, transcript_id, transcript_type) %>%
    dplyr::distinct()
  
  annotated_orfs <- annotate_orf_isoforms(
    annotations, 
    orfs,
    BSgenome,
    transcripts_meta
  )
  
  # Create mock DTU results
  test_transcripts <- unique(annotated_orfs$table$transcript_id)[1:20]
  test_genes <- unique(annotated_orfs$table$gene_id[
    annotated_orfs$table$transcript_id %in% test_transcripts
  ])
  
  dtu_results <- data.frame(
    transcript_id = test_transcripts,
    gene_id = rep(test_genes, length.out = length(test_transcripts)),
    gene_name = paste0("GENE", seq_along(test_transcripts)),
    pvalue = runif(length(test_transcripts)),
    adj_pvalue = runif(length(test_transcripts)),
    rank = rnorm(length(test_transcripts))
  )
  
  selected_genes <- test_genes[1:min(3, length(test_genes))]
  
  # Test function execution
  expect_no_error({
    results <- diff_orf_usage(
      annotated_orfs$table,
      dtu_results,
      selected_genes
    )
  })
  
  # Run function
  results <- diff_orf_usage(
    annotated_orfs$table,
    dtu_results,
    selected_genes
  )
  
  # Test output structure
  expect_type(results, "list")
  expect_named(results, c("gene_stats", "gene_tables"))
  
  # Test gene_stats structure
  expect_s3_class(results$gene_stats, "data.frame")
  expect_true(nrow(results$gene_stats) > 0)
  
  required_stat_cols <- c("gene_id", "gene_name", "n_tx_isoforms", 
                          "n_signif_dtu", "n_orfs", "n_events", "n_translatable")
  expect_true(all(required_stat_cols %in% colnames(results$gene_stats)))
  
  # Test gene_tables structure
  expect_type(results$gene_tables, "list")
  expect_equal(length(results$gene_tables), nrow(results$gene_stats))
  expect_true(all(names(results$gene_tables) %in% results$gene_stats$gene_id))
  
  # Test individual gene table structure
  first_gene_table <- results$gene_tables[[1]]
  expect_s3_class(first_gene_table, "data.table")
  expect_true("transcript_id" %in% colnames(first_gene_table))
  expect_true("adj_pvalue" %in% colnames(first_gene_table))
  expect_true("rank" %in% colnames(first_gene_table))
  
  # Test that numeric columns are valid
  expect_true(all(results$gene_stats$n_tx_isoforms > 0))
  expect_true(all(results$gene_stats$n_signif_dtu >= 0))
  expect_true(all(results$gene_stats$n_orfs > 0))
  expect_true(all(results$gene_stats$n_events >= 0))
  expect_true(all(results$gene_stats$n_translatable >= 0))
  
  # Test sorting by n_events
  if (nrow(results$gene_stats) > 1) {
    expect_true(all(diff(results$gene_stats$n_events) <= 0))
  }
})

test_that("diff_orf_usage validates required columns", {
  # Test with missing required columns
  incomplete_dtu <- data.frame(
    transcript_id = c("tx1", "tx2"),
    gene_id = c("gene1", "gene1")
    # Missing: gene_name, pvalue, rank, adj_pvalue
  )
  
  mock_orf_tab <- data.frame(
    ORF_id = c("orf1", "orf2"),
    transcript_id = c("tx1", "tx2"),
    gene_id = c("gene1", "gene1"),
    unique_tx_iso = c("translatable_1", "translatable_1")
  )
  
  expect_error(
    diff_orf_usage(mock_orf_tab, incomplete_dtu, "gene1"),
    "Not all required columns are present"
  )
})

test_that("diff_orf_usage handles genes with no ORFs", {
  # Create DTU results for genes without ORFs
  dtu_results <- data.frame(
    transcript_id = c("tx1", "tx2", "tx3"),
    gene_id = c("gene1", "gene1", "gene2"),
    gene_name = c("GENEA", "GENEA", "GENEB"),
    pvalue = c(0.01, 0.05, 0.001),
    adj_pvalue = c(0.02, 0.1, 0.005),
    rank = c(-1.5, 2.0, -2.5)
  )
  
  # ORF table with different genes
  orf_tab <- data.frame(
    ORF_id = c("orf1", "orf2"),
    transcript_id = c("tx4", "tx5"),
    gene_id = c("gene3", "gene3"),
    unique_tx_iso = c("translatable_1", "translatable_1")
  )
  
  # This should run but return empty or minimal results
  results <- diff_orf_usage(orf_tab, dtu_results, c("gene1", "gene2"))
  
  expect_type(results, "list")
  expect_named(results, c("gene_stats", "gene_tables"))
})

test_that("diff_orf_usage handles single gene analysis", {
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
  skip_if(gtf == "" || bed == "", "Test data not found")
  
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(rtracklayer))
  BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  
  annotations <- suppressWarnings(suppressMessages(
    prepare_annotations_fromGTF(gtf, BSgenome)
  ))
  
  orfs <- import(bed, format = "BED")
  names(orfs) <- orfs$name
  
  transcripts_meta <- import(gtf, format = "GTF") %>%
    as.data.frame() %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::select(gene_name, gene_id, transcript_id, transcript_type) %>%
    dplyr::distinct()
  
  annotated_orfs <- annotate_orf_isoforms(
    annotations, 
    orfs[1:10],
    BSgenome,
    transcripts_meta
  )
  
  # Get one gene
  test_gene <- unique(annotated_orfs$table$gene_id)[1]
  test_transcripts <- annotated_orfs$table$transcript_id[
    annotated_orfs$table$gene_id == test_gene
  ]
  
  dtu_results <- data.frame(
    transcript_id = test_transcripts,
    gene_id = test_gene,
    gene_name = "TEST_GENE",
    pvalue = runif(length(test_transcripts)),
    adj_pvalue = runif(length(test_transcripts)),
    rank = rnorm(length(test_transcripts))
  )
  
  # Analyze single gene
  results <- diff_orf_usage(
    annotated_orfs$table,
    dtu_results,
    test_gene
  )
  
  # Should have exactly one gene in results
  expect_equal(nrow(results$gene_stats), 1)
  expect_equal(length(results$gene_tables), 1)
  expect_equal(results$gene_stats$gene_id[1], test_gene)
})

test_that("diff_orf_usage correctly identifies translatable events", {
  # Create a simple test case
  orf_tab <- data.frame(
    ORF_id = c("orf1", "orf1", "orf2", "orf2"),
    transcript_id = c("tx1", "tx2", "tx1", "tx2"),
    gene_id = c("gene1", "gene1", "gene1", "gene1"),
    unique_tx_iso = c("translatable_1", "internal_stop", 
                      "translatable_1", "translatable_2")
  )
  
  dtu_results <- data.frame(
    transcript_id = c("tx1", "tx2"),
    gene_id = c("gene1", "gene1"),
    gene_name = c("GENEA", "GENEA"),
    pvalue = c(0.01, 0.05),
    adj_pvalue = c(0.02, 0.1),
    rank = c(-1.5, 2.0)
  )
  
  results <- diff_orf_usage(orf_tab, dtu_results, "gene1")
  
  # Check that translatable events are counted
  expect_true(results$gene_stats$n_translatable > 0)
  expect_true(results$gene_stats$n_events >= results$gene_stats$n_translatable)
})

test_that("diff_orf_usage counts significant DTU events", {
  orf_tab <- data.frame(
    ORF_id = rep("orf1", 3),
    transcript_id = c("tx1", "tx2", "tx3"),
    gene_id = rep("gene1", 3),
    unique_tx_iso = rep("translatable_1", 3)
  )
  
  dtu_results <- data.frame(
    transcript_id = c("tx1", "tx2", "tx3"),
    gene_id = rep("gene1", 3),
    gene_name = rep("GENEA", 3),
    pvalue = c(0.001, 0.01, 0.5),
    adj_pvalue = c(0.01, 0.04, 0.6),  # 2 significant at 0.05 level
    rank = c(-2, -1, 0.5)
  )
  
  results <- diff_orf_usage(orf_tab, dtu_results, "gene1")
  
  # Should count 2 significant DTU events (adj_pvalue < 0.05)
  expect_equal(results$gene_stats$n_signif_dtu, 2)
})