# Helper function to create test GRanges
create_test_gr <- function(coords, tx_ids = NULL, orf_types = NULL, orf_ids = NULL) {
  suppressMessages(library(GenomicRanges))
  
  gr <- GRanges(
    seqnames = coords$chr,
    ranges = IRanges(start = coords$start, end = coords$end),
    strand = coords$strand
  )
  
  if (!is.null(tx_ids)) {
    mcols(gr)$transcript_id <- tx_ids
  }
  if (!is.null(orf_types)) {
    mcols(gr)$orf_type <- orf_types
  }
  if (!is.null(orf_ids)) {
    mcols(gr)$orf_id <- orf_ids
    names(gr) <- orf_ids
  }
  
  gr
}

# =============================================================================
# Basic functionality tests
# =============================================================================

test_that("collapse_orf_calls works with exact coordinate matching", {
  suppressMessages(library(GenomicRanges))
  
  # Create three datasets with some overlapping ORFs
  gr1 <- create_test_gr(
    data.frame(chr = c("chr1", "chr1", "chr2"),
               start = c(100, 500, 200),
               end = c(300, 800, 400),
               strand = c("+", "-", "+")),
    tx_ids = c("tx1", "tx2", "tx3"),
    orf_types = c("uORF", "CDS", "dORF"),
    orf_ids = c("orf1", "orf2", "orf3")
  )
  
  gr2 <- create_test_gr(
    data.frame(chr = c("chr1", "chr2"),
               start = c(100, 200),
               end = c(300, 400),
               strand = c("+", "+")),
    tx_ids = c("tx1", "tx3"),
    orf_types = c("uORF", "dORF"),
    orf_ids = c("orf1", "orf3")
  )
  
  gr3 <- create_test_gr(
    data.frame(chr = c("chr1"),
               start = c(500),
               end = c(800),
               strand = c("-")),
    tx_ids = c("tx2"),
    orf_types = c("CDS"),
    orf_ids = c("orf2")
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2, sample3 = gr3)
  
  result <- collapse_orf_calls(orf_list, match_by = "coordinates")
  
  # Should have 3 unique ORFs (union of all)
  expect_equal(length(result), 3)
  
  # Check presence columns exist
  expect_true("in_sample1" %in% colnames(mcols(result)))
  expect_true("in_sample2" %in% colnames(mcols(result)))
  expect_true("in_sample3" %in% colnames(mcols(result)))
  expect_true("n_datasets" %in% colnames(mcols(result)))
  
  # Check ORF1 (chr1:100-300:+) is in samples 1 and 2
  orf1_idx <- which(seqnames(result) == "chr1" & start(result) == 100)
  expect_true(mcols(result)$in_sample1[orf1_idx])
  expect_true(mcols(result)$in_sample2[orf1_idx])
  expect_false(mcols(result)$in_sample3[orf1_idx])
  expect_equal(mcols(result)$n_datasets[orf1_idx], 2)
  
  # Check ORF2 (chr1:500-800:-) is in samples 1 and 3
  orf2_idx <- which(seqnames(result) == "chr1" & start(result) == 500)
  expect_true(mcols(result)$in_sample1[orf2_idx])
  expect_false(mcols(result)$in_sample2[orf2_idx])
  expect_true(mcols(result)$in_sample3[orf2_idx])
  expect_equal(mcols(result)$n_datasets[orf2_idx], 2)
  
  # Check ORF3 (chr2:200-400:+) is in samples 1 and 2
  orf3_idx <- which(seqnames(result) == "chr2")
  expect_true(mcols(result)$in_sample1[orf3_idx])
  expect_true(mcols(result)$in_sample2[orf3_idx])
  expect_false(mcols(result)$in_sample3[orf3_idx])
  expect_equal(mcols(result)$n_datasets[orf3_idx], 2)
})

test_that("collapse_orf_calls drops ORF IDs by default", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    orf_ids = "orf1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    orf_ids = "orf2"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list)
  
  # Names should be NULL
  expect_null(names(result))
  
  # orf_id column should not exist
  expect_false("orf_id" %in% colnames(mcols(result)))
})

test_that("collapse_orf_calls keeps ORF IDs when requested", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    orf_ids = "orf1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 500, end = 800, strand = "+"),
    orf_ids = "orf2"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, keep_orf_ids = TRUE)
  
  # Names should be present
  expect_false(is.null(names(result)))
  expect_equal(length(names(result)), 2)
  
  # orf_id column should exist
  expect_true("orf_id" %in% colnames(mcols(result)))
})

# =============================================================================
# Match strategy tests
# =============================================================================

test_that("match_by = 'start_codon' matches by start position only", {
  suppressMessages(library(GenomicRanges))
  
  # Same start, different stops
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 400, strand = "+"),
    tx_ids = "tx2"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, match_by = "start_codon")
  
  # Should collapse to 1 ORF (same start)
  expect_equal(length(result), 1)
  expect_equal(mcols(result)$n_datasets[1], 2)
})

test_that("match_by = 'start_codon_transcript' requires same transcript", {
  suppressMessages(library(GenomicRanges))
  
  # Same start, different transcripts
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 400, strand = "+"),
    tx_ids = "tx2"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, match_by = "start_codon_transcript")
  
  # Should NOT collapse (different transcripts)
  expect_equal(length(result), 2)
  expect_equal(mcols(result)$n_datasets[1], 1)
  expect_equal(mcols(result)$n_datasets[2], 1)
})

test_that("match_by = 'stop_codon' matches by stop position only", {
  suppressMessages(library(GenomicRanges))
  
  # Different starts, same stop
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 50, end = 300, strand = "+"),
    tx_ids = "tx2"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, match_by = "stop_codon")
  
  # Should collapse to 1 ORF (same stop)
  expect_equal(length(result), 1)
  expect_equal(mcols(result)$n_datasets[1], 2)
})

test_that("match_by = 'stop_codon_transcript' requires same transcript", {
  suppressMessages(library(GenomicRanges))
  
  # Same stop, different transcripts
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 50, end = 300, strand = "+"),
    tx_ids = "tx2"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, match_by = "stop_codon_transcript")
  
  # Should NOT collapse (different transcripts)
  expect_equal(length(result), 2)
})

test_that("match_by = 'start_stop' ignores transcript differences", {
  suppressMessages(library(GenomicRanges))
  
  # Same coordinates, different transcripts
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx2"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, match_by = "start_stop")
  
  # Should collapse to 1 ORF
  expect_equal(length(result), 1)
  expect_equal(mcols(result)$n_datasets[1], 2)
})

test_that("match_by = 'start_stop_transcript' requires exact match", {
  suppressMessages(library(GenomicRanges))
  
  # Same coordinates, different transcripts
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx2"
  )
  gr3 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2, sample3 = gr3)
  
  result <- collapse_orf_calls(orf_list, match_by = "start_stop_transcript")
  
  # Should have 2 ORFs (tx1 and tx2)
  expect_equal(length(result), 2)
  
  # tx1 should be in samples 1 and 3
  tx1_orf <- result[mcols(result)$transcript_id == "tx1"]
  expect_equal(mcols(tx1_orf)$n_datasets[1], 2)
  expect_true(mcols(tx1_orf)$in_sample1[1])
  expect_false(mcols(tx1_orf)$in_sample2[1])
  expect_true(mcols(tx1_orf)$in_sample3[1])
})

test_that("match_by = 'overlap' matches overlapping ORFs", {
  suppressMessages(library(GenomicRanges))
  
  # Slightly different coordinates with good overlap
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 105, end = 305, strand = "+")
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, match_by = "overlap", 
                               overlap_threshold = 0.9)
  
  # Should collapse to 1 ORF (high overlap)
  expect_equal(length(result), 1)
  expect_equal(mcols(result)$n_datasets[1], 2)
})

test_that("match_by = 'metadata' matches by custom keys", {
  suppressMessages(library(GenomicRanges))
  
  # Different coordinates but same metadata
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1",
    orf_types = "uORF"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 500, end = 700, strand = "+"),
    tx_ids = "tx1",
    orf_types = "uORF"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, match_by = "metadata",
                               metadata_keys = c("transcript_id", "orf_type"))
  
  # Should collapse to 1 ORF (same transcript + orf_type)
  expect_equal(length(result), 1)
  expect_equal(mcols(result)$n_datasets[1], 2)
})

# =============================================================================
# Combine strategy tests
# =============================================================================

test_that("combine = 'first' uses first dataset as reference", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = c("chr1", "chr2"),
               start = c(100, 200),
               end = c(300, 400),
               strand = c("+", "+"))
  )
  gr2 <- create_test_gr(
    data.frame(chr = c("chr1", "chr3"),
               start = c(100, 500),
               end = c(300, 700),
               strand = c("+", "+"))
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, combine = "first")
  
  # Should only have ORFs from first dataset
  expect_equal(length(result), 2)
  
  # All should be in sample1
  expect_true(all(mcols(result)$in_sample1))
  
  # chr3 ORF should not be present (only in sample2)
  expect_false(any(as.character(seqnames(result)) == "chr3"))
})

test_that("combine = 'union' includes all ORFs", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = c("chr1"),
               start = c(100),
               end = c(300),
               strand = c("+")),
    tx_ids = "tx1"
  )
  gr2 <- create_test_gr(
    data.frame(chr = c("chr2"),
               start = c(500),
               end = c(700),
               strand = c("+")),
    tx_ids = "tx2"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list, combine = "union")
  
  # Should have ORFs from both datasets
  expect_equal(length(result), 2)
  expect_true("chr1" %in% as.character(seqnames(result)))
  expect_true("chr2" %in% as.character(seqnames(result)))
})

# =============================================================================
# Edge cases and error handling
# =============================================================================

test_that("collapse_orf_calls handles single dataset", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  
  expect_warning(
    result <- collapse_orf_calls(list(sample1 = gr1)),
    "Only one dataset provided"
  )
  
  expect_equal(length(result), 1)
  expect_true("n_datasets" %in% colnames(mcols(result)))
  expect_equal(mcols(result)$n_datasets[1], 1)
  expect_true("in_sample1" %in% colnames(mcols(result)))
})

test_that("collapse_orf_calls errors with empty list", {
  expect_error(
    collapse_orf_calls(list()),
    "orf_list must be a list with at least one element"
  )
})

test_that("collapse_orf_calls errors with non-GRanges input", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  
  expect_error(
    collapse_orf_calls(list(gr1, data.frame(a = 1))),
    "All elements of orf_list must be GRanges objects"
  )
})

test_that("collapse_orf_calls uses list names for dataset names", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  
  orf_list <- list(control = gr1, treatment = gr2)
  
  result <- collapse_orf_calls(orf_list)
  
  expect_true("in_control" %in% colnames(mcols(result)))
  expect_true("in_treatment" %in% colnames(mcols(result)))
})

test_that("collapse_orf_calls generates default dataset names", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  
  orf_list <- list(gr1, gr2)  # No names
  
  result <- collapse_orf_calls(orf_list)
  
  expect_true("in_dataset_1" %in% colnames(mcols(result)))
  expect_true("in_dataset_2" %in% colnames(mcols(result)))
})

test_that("collapse_orf_calls handles custom dataset names", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  
  orf_list <- list(gr1, gr2)
  
  result <- collapse_orf_calls(orf_list, dataset_names = c("replicate_1", "replicate_2"))
  
  expect_true("in_replicate_1" %in% colnames(mcols(result)))
  expect_true("in_replicate_2" %in% colnames(mcols(result)))
})

test_that("collapse_orf_calls errors when transcript_id missing for transcript-based matching", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+")
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  expect_error(
    collapse_orf_calls(orf_list, match_by = "start_codon_transcript"),
    "transcript_id not found in all datasets"
  )
})

# =============================================================================
# Metadata preservation tests
# =============================================================================

test_that("collapse_orf_calls preserves metadata columns", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = "chr1", start = 100, end = 300, strand = "+"),
    tx_ids = "tx1",
    orf_types = "uORF"
  )
  gr2 <- create_test_gr(
    data.frame(chr = "chr1", start = 500, end = 700, strand = "+"),
    tx_ids = "tx2",
    orf_types = "CDS"
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2)
  
  result <- collapse_orf_calls(orf_list)
  
  # Metadata columns should be preserved
  expect_true("transcript_id" %in% colnames(mcols(result)))
  expect_true("orf_type" %in% colnames(mcols(result)))
})

test_that("n_datasets count is accurate", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- create_test_gr(
    data.frame(chr = c("chr1", "chr2", "chr3"),
               start = c(100, 200, 300),
               end = c(150, 250, 350),
               strand = c("+", "+", "+"))
  )
  gr2 <- create_test_gr(
    data.frame(chr = c("chr1", "chr2"),
               start = c(100, 200),
               end = c(150, 250),
               strand = c("+", "+"))
  )
  gr3 <- create_test_gr(
    data.frame(chr = c("chr1"),
               start = c(100),
               end = c(150),
               strand = c("+"))
  )
  
  orf_list <- list(sample1 = gr1, sample2 = gr2, sample3 = gr3)
  
  result <- collapse_orf_calls(orf_list)
  
  # chr1 should be in all 3
  chr1_orf <- result[seqnames(result) == "chr1" & start(result) == 100]
  expect_equal(mcols(chr1_orf)$n_datasets[1], 3)
  
  # chr2 should be in 2
  chr2_orf <- result[seqnames(result) == "chr2"]
  expect_equal(mcols(chr2_orf)$n_datasets[1], 2)
  
  # chr3 should be in 1
  chr3_orf <- result[seqnames(result) == "chr3"]
  expect_equal(mcols(chr3_orf)$n_datasets[1], 1)
})
