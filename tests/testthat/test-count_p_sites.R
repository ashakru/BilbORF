test_that("count_p_sites quantifies P-sites correctly", {
  # Setup: Load test data
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
  bam <- system.file("extdata", "riboseq_chr10.bam", package = "BilbORF")
  skip_if(gtf == "" || bed == "" || bam == "", "Test data not found")
  
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(rtracklayer))
  BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  
  # Prepare annotations and ORFs
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
    orfs[1:20],  # Use subset for faster testing
    BSgenome,
    transcripts_meta
  )
  
  # Define offsets
  offsets <- data.frame(
    fraction = 25:30,
    offsets_start = -12
  )
  
  # Test function execution
  expect_no_error({
    results <- count_p_sites(bam, offsets, annotated_orfs)
  })
  
  # Run function
  results <- count_p_sites(bam, offsets, annotated_orfs)
  
  # Test output structure
  expect_s3_class(results, "data.frame")
  expect_true(nrow(results) > 0)
  
  # Test required columns
  expected_cols <- c("ORF_id", "entropy", "frame_0", "frame_1", "frame_2", "frames_sum")
  expect_true(all(expected_cols %in% colnames(results)))
  
  # Test that coverage values are non-negative
  expect_true(all(results$frame_0 >= 0, na.rm = TRUE))
  expect_true(all(results$frame_1 >= 0, na.rm = TRUE))
  expect_true(all(results$frame_2 >= 0, na.rm = TRUE))
  
  # Test that frames_sum equals sum of individual frames
  calculated_sum <- results$frame_0 + results$frame_1 + results$frame_2
  expect_equal(results$frames_sum, calculated_sum)
  
  # Test entropy values are valid (between 0 and log2(3) for 3 frames)
  expect_true(all(results$entropy >= 0, na.rm = TRUE))
  expect_true(all(results$entropy <= log2(3), na.rm = TRUE))
})

test_that("count_p_sites handles different offset configurations", {
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
  bam <- system.file("extdata", "riboseq_chr10.bam", package = "BilbORF")
  skip_if(gtf == "" || bed == "" || bam == "", "Test data not found")
  
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
  
  # Test with single read length
  offsets1 <- data.frame(fraction = 28, offsets_start = -12)
  expect_no_error({
    results1 <- count_p_sites(bam, offsets1, annotated_orfs)
  })
  
  # Test with multiple read lengths
  offsets2 <- data.frame(
    fraction = 26:30,
    offsets_start = -12
  )
  expect_no_error({
    results2 <- count_p_sites(bam, offsets2, annotated_orfs)
  })
  
  # Both should return valid data frames
  expect_s3_class(results1, "data.frame")
  expect_s3_class(results2, "data.frame")
})

test_that("count_p_sites validates inputs", {
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
    orfs[1:5],
    BSgenome,
    transcripts_meta
  )
  
  offsets <- data.frame(fraction = 25:30, offsets_start = -12)
  
  # Test with non-existent BAM file
  expect_error(
    count_p_sites("nonexistent.bam", offsets, annotated_orfs)
  )
  
  # Test with NULL offsets
  bam <- system.file("extdata", "riboseq_chr10.bam", package = "BilbORF")
  skip_if(bam == "", "BAM file not found")
  
  expect_error(
    count_p_sites(bam, NULL, annotated_orfs)
  )
  
  # Test with NULL annotated_orfs
  expect_error(
    count_p_sites(bam, offsets, NULL)
  )
})

test_that("count_p_sites handles ORFs with missing coverage", {
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
  bam <- system.file("extdata", "riboseq_chr10.bam", package = "BilbORF")
  skip_if(gtf == "" || bed == "" || bam == "", "Test data not found")
  
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
    orfs[1:30],
    BSgenome,
    transcripts_meta
  )
  
  offsets <- data.frame(fraction = 25:30, offsets_start = -12)
  results <- count_p_sites(bam, offsets, annotated_orfs)
  
  # Should handle ORFs with zero coverage gracefully
  # Check if any ORFs have zero total coverage
  zero_coverage <- results$frames_sum == 0
  if (any(zero_coverage)) {
    # These should still have valid (possibly 0) values for all metrics
    expect_true(all(!is.na(results$frame_0[zero_coverage])))
    expect_true(all(!is.na(results$frame_1[zero_coverage])))
    expect_true(all(!is.na(results$frame_2[zero_coverage])))
  }
})