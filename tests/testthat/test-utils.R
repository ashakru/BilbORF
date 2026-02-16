test_that("rank_exons assigns correct ranks for positive strand", {
  suppressMessages(library(GenomicRanges))
  
  # Create test GRangesList with positive strand
  gr_pos <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(100, 200, 300), end = c(150, 250, 350)),
    strand = "+"
  )
  grl <- GRangesList(tx1 = gr_pos)
  
  # Apply rank_exons
  ranked <- rank_exons(grl)
  
  # Check that ranks go from 1 to 3 for positive strand
  expect_equal(ranked$tx1$exon_rank, c(1, 2, 3))
  
  # Check that ranges are sorted by rank
  expect_true(all(start(ranked$tx1) == c(100, 200, 300)))
})

test_that("rank_exons assigns correct ranks for negative strand", {
  suppressMessages(library(GenomicRanges))
  
  # Create test GRangesList with negative strand
  gr_neg <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(100, 200, 300), end = c(150, 250, 350)),
    strand = "-"
  )
  grl <- GRangesList(tx1 = gr_neg)
  
  # Apply rank_exons
  ranked <- rank_exons(grl)
  
  # Check that ranks go from 3 to 1 for negative strand (reverse order)
  expect_equal(ranked$tx1$exon_rank, c(3, 2, 1))
  
  # After sorting by rank, should be in reverse order
  expect_true(all(start(ranked$tx1) == c(300, 200, 100)))
})

test_that("rank_exons handles multiple transcripts", {
  suppressMessages(library(GenomicRanges))
  
  gr_pos <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(100, 200), end = c(150, 250)),
    strand = "+"
  )
  
  gr_neg <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(300, 400, 500), end = c(350, 450, 550)),
    strand = "-"
  )
  
  grl <- GRangesList(tx_pos = gr_pos, tx_neg = gr_neg)
  
  # Apply rank_exons
  ranked <- rank_exons(grl)
  
  # Check positive strand transcript
  expect_equal(ranked$tx_pos$exon_rank, c(1, 2))
  
  # Check negative strand transcript
  expect_equal(ranked$tx_neg$exon_rank, c(3, 2, 1))
})

test_that("rank_exons returns GRangesList", {
  suppressMessages(library(GenomicRanges))
  
  gr <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(100, 200), end = c(150, 250)),
    strand = "+"
  )
  grl <- GRangesList(tx1 = gr)
  
  ranked <- rank_exons(grl)
  
  expect_s4_class(ranked, "GRangesList")
  expect_equal(length(ranked), 1)
  expect_true("exon_rank" %in% colnames(mcols(ranked$tx1)))
})

test_that("get_start_position returns minimum start coordinate", {
  suppressMessages(library(GenomicRanges))
  
  # Create GRangesList with multiple ranges
  gr1 <- GRanges("chr1", IRanges(c(100, 200, 300), c(150, 250, 350)), strand = "+")
  gr2 <- GRanges("chr1", IRanges(c(500, 600), c(550, 650)), strand = "+")
  grl <- GRangesList(orf1 = gr1, orf2 = gr2)
  
  starts <- get_start_position(grl)
  
  expect_equal(starts["orf1"], 100)
  expect_equal(starts["orf2"], 500)
  expect_length(starts, 2)
})

test_that("get_stop_position returns maximum end coordinate", {
  suppressMessages(library(GenomicRanges))
  
  # Create GRangesList with multiple ranges
  gr1 <- GRanges("chr1", IRanges(c(100, 200, 300), c(150, 250, 350)), strand = "+")
  gr2 <- GRanges("chr1", IRanges(c(500, 600), c(550, 650)), strand = "+")
  grl <- GRangesList(orf1 = gr1, orf2 = gr2)
  
  ends <- get_stop_position(grl)
  
  expect_equal(ends["orf1"], 350)
  expect_equal(ends["orf2"], 650)
  expect_length(ends, 2)
})

test_that("get_start_position and get_stop_position handle single exon ORFs", {
  suppressMessages(library(GenomicRanges))
  
  gr <- GRanges("chr1", IRanges(100, 200), strand = "+")
  grl <- GRangesList(orf1 = gr)
  
  starts <- get_start_position(grl)
  ends <- get_stop_position(grl)
  
  expect_equal(starts["orf1"], 100)
  expect_equal(ends["orf1"], 200)
})

test_that("get_start_position handles empty GRangesList", {
  suppressMessages(library(GenomicRanges))
  
  grl <- GRangesList()
  starts <- get_start_position(grl)
  
  expect_length(starts, 0)
  expect_type(starts, "integer")
})

test_that("utility functions preserve names", {
  suppressMessages(library(GenomicRanges))
  
  gr1 <- GRanges("chr1", IRanges(c(100, 200), c(150, 250)), strand = "+")
  gr2 <- GRanges("chr2", IRanges(c(300, 400), c(350, 450)), strand = "-")
  grl <- GRangesList(gene1_orf1 = gr1, gene2_orf2 = gr2)
  
  starts <- get_start_position(grl)
  ends <- get_stop_position(grl)
  
  expect_equal(names(starts), c("gene1_orf1", "gene2_orf2"))
  expect_equal(names(ends), c("gene1_orf1", "gene2_orf2"))
})

test_that("rank_exons handles mixed strand GRangesList", {
  suppressMessages(library(GenomicRanges))
  
  # Should handle case where some exons might have * strand
  gr_unstranded <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(100, 200), end = c(150, 250)),
    strand = "*"
  )
  grl <- GRangesList(tx1 = gr_unstranded)
  
  # Should still assign ranks (treated as positive)
  ranked <- rank_exons(grl)
  expect_true("exon_rank" %in% colnames(mcols(ranked$tx1)))
  expect_equal(ranked$tx1$exon_rank, c(1, 2))
})

test_that("rank_exons maintains other metadata", {
  suppressMessages(library(GenomicRanges))
  
  gr <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(100, 200), end = c(150, 250)),
    strand = "+",
    gene_id = c("ENSG001", "ENSG001"),
    transcript_id = c("ENST001", "ENST001")
  )
  grl <- GRangesList(tx1 = gr)
  
  ranked <- rank_exons(grl)
  
  # Original metadata should be preserved
  expect_true("gene_id" %in% colnames(mcols(ranked$tx1)))
  expect_true("transcript_id" %in% colnames(mcols(ranked$tx1)))
  expect_equal(ranked$tx1$gene_id, c("ENSG001", "ENSG001"))
})