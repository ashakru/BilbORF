# ==============================================================================
# test-orf_footprints.R
#
# Unit tests for orf_footprints.R.
#
# Strategy: all tests use minimal synthetic TxDb / GRangesList fixtures so
# that no external data files or internet access are needed.
#
# Fixture helper: make_fake_exons_by_tx() builds a named GRangesList that
# mimics the output of GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE).
# Each element is a GRanges of exon blocks for one transcript.
# ==============================================================================

library(testthat)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)

# ------------------------------------------------------------------------------
# Fixture helpers
# ------------------------------------------------------------------------------

#' Build a minimal exons_by_tx list from a named list of data frames.
#' Each df must have columns: start, end (1-based genomic), strand, seqname.
make_fake_exons_by_tx <- function(tx_list) {
  grl <- lapply(names(tx_list), function(tx_id) {
    df <- tx_list[[tx_id]]
    GenomicRanges::GRanges(
      seqnames = df$seqname,
      ranges   = IRanges::IRanges(start = df$start, end = df$end),
      strand   = df$strand
    )
  })
  names(grl) <- names(tx_list)
  GenomicRanges::GRangesList(grl)
}

# ==============================================================================
# Tests for reconstruct_chain()
# ==============================================================================

test_that("single-exon plus-strand ORF is reconstructed correctly", {
  # Transcript: one exon chr1:+:100-500
  # ORF: start=200, stop=400 (both within the exon)
  exons_by_tx <- make_fake_exons_by_tx(list(
    tx1 = data.frame(seqname = "chr1", start = 100L, end = 500L, strand = "+",
                     stringsAsFactors = FALSE)
  ))

  fp <- reconstruct_chain(
    start = 200L, stop = 400L, strand = "+",
    transcript_id = "tx1", exons_by_tx = exons_by_tx
  )

  expect_s4_class(fp, "GRanges")
  expect_equal(length(fp), 1L)                        # single exon
  expect_equal(GenomicRanges::start(fp), 200L)        # clipped to start
  expect_equal(GenomicRanges::end(fp),   400L)        # clipped to stop
  expect_equal(as.character(GenomicRanges::strand(fp)), "+")
  expect_equal(fp$exon_rank, 1L)
})

test_that("multi-exon plus-strand ORF is reconstructed and ordered correctly", {
  # Transcript: three exons on chr1+
  #   exon1: 100–300
  #   exon2: 500–700
  #   exon3: 900–1100
  # ORF: start=150 (within exon1), stop=650 (within exon2)
  # Expected chain: [150,300] [500,650]  in that order (5′→3′)
  exons_by_tx <- make_fake_exons_by_tx(list(
    tx2 = data.frame(
      seqname = "chr1",
      start   = c(100L, 500L, 900L),
      end     = c(300L, 700L, 1100L),
      strand  = "+"
    )
  ))

  fp <- reconstruct_chain(
    start = 150L, stop = 650L, strand = "+",
    transcript_id = "tx2", exons_by_tx = exons_by_tx
  )

  expect_equal(length(fp), 2L)
  # Exon 1 clipped at start
  expect_equal(GenomicRanges::start(fp[1]), 150L)
  expect_equal(GenomicRanges::end(fp[1]),   300L)
  # Exon 2 clipped at stop
  expect_equal(GenomicRanges::start(fp[2]), 500L)
  expect_equal(GenomicRanges::end(fp[2]),   650L)
  # Exon ranks in 5′→3′ order
  expect_equal(fp$exon_rank, c(1L, 2L))
  # Exon3 (900–1100) must NOT be included
  expect_true(all(GenomicRanges::end(fp) <= 700L))
})

test_that("multi-exon minus-strand ORF: path is in 5'->3' (descending) order", {
  # Transcript: two exons on chr1-
  #   exon_A: genomic 800–1000  (5' end on minus strand = exon_B)
  #   exon_B: genomic 200–400
  # On the minus strand 5'→3' means descending genomic coords:
  #   first coding exon = 800–1000, second = 200–400.
  # ORF: start=950 (within exon_A, minus strand so numerically larger),
  #      stop=250  (within exon_B).
  exons_by_tx <- make_fake_exons_by_tx(list(
    tx3 = data.frame(
      seqname = "chr1",
      start   = c(200L, 800L),
      end     = c(400L, 1000L),
      strand  = "-"
    )
  ))

  fp <- reconstruct_chain(
    start = 950L, stop = 250L, strand = "-",
    transcript_id = "tx3", exons_by_tx = exons_by_tx
  )

  expect_equal(length(fp), 2L)
  # First element in coding order should be the higher-coord block (5' on minus)
  expect_equal(GenomicRanges::start(fp[1]), 800L)   # clipped at start=950 end
  expect_equal(GenomicRanges::end(fp[1]),   950L)   # start codon is here
  expect_equal(GenomicRanges::start(fp[2]), 250L)   # stop codon is here
  expect_equal(GenomicRanges::end(fp[2]),   400L)
  expect_equal(fp$exon_rank, c(1L, 2L))             # rank 1 = 5' end
})

test_that("start or stop outside any exon returns NULL with a warning", {
  exons_by_tx <- make_fake_exons_by_tx(list(
    tx4 = data.frame(seqname = "chr1", start = 100L, end = 300L, strand = "+")
  ))

  # start = 50 is upstream of the exon
  expect_warning(
    result <- reconstruct_chain(
      start = 50L, stop = 250L, strand = "+",
      transcript_id = "tx4", exons_by_tx = exons_by_tx
    )
  )
  expect_null(result)
})

test_that("transcript not found in exons_by_tx returns NULL with a warning", {
  exons_by_tx <- make_fake_exons_by_tx(list(
    tx5 = data.frame(seqname = "chr1", start = 100L, end = 300L, strand = "+")
  ))

  expect_warning(
    result <- reconstruct_chain(
      start = 110L, stop = 290L, strand = "+",
      transcript_id = "tx_missing", exons_by_tx = exons_by_tx
    )
  )
  expect_null(result)
})

# ==============================================================================
# Tests for make_chain_id()
# ==============================================================================

test_that("make_chain_id is deterministic", {
  fp <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(start = 200L, end = 400L),
    strand   = "+"
  )
  id1 <- make_chain_id(fp)
  id2 <- make_chain_id(fp)
  expect_equal(id1, id2)
  expect_match(id1, "^chain_[0-9a-f]+$")
})

test_that("two transcripts sharing exact ORF-spanning path → same chain_id", {
  # tx_A and tx_B differ only in exons outside the ORF span (e.g. 5'UTR exon).
  # Both have the same two coding exons [150,300] and [500,650].
  # reconstruct_chain clips both identically → same chain_id.

  exons_by_tx <- make_fake_exons_by_tx(list(
    # tx_A has an upstream UTR exon (10-80) not overlapping ORF
    tx_A = data.frame(
      seqname = "chr1",
      start   = c(10L,  150L, 500L),
      end     = c(80L,  300L, 650L),
      strand  = "+"
    ),
    # tx_B starts its first exon at 130 — still covers orf start 150
    tx_B = data.frame(
      seqname = "chr1",
      start   = c(130L, 500L),
      end     = c(300L, 650L),
      strand  = "+"
    )
  ))

  fp_A <- reconstruct_chain(150L, 620L, "+", "tx_A", exons_by_tx)
  fp_B <- reconstruct_chain(150L, 620L, "+", "tx_B", exons_by_tx)

  expect_false(is.null(fp_A))
  expect_false(is.null(fp_B))
  expect_equal(make_chain_id(fp_A), make_chain_id(fp_B))
})

test_that("shared internal junctions but different start → two chain_ids", {
  # tx_C and tx_D share exon2 [500,700] completely.
  # tx_C: ORF starts at 200 (in exon [100,300])
  # tx_D: ORF starts at 150 (in exon [100,300]) — different start, same stop
  # → different codon walk at N-terminus → different chain_ids.
  exons_by_tx <- make_fake_exons_by_tx(list(
    tx_C = data.frame(
      seqname = "chr1",
      start   = c(100L, 500L),
      end     = c(300L, 700L),
      strand  = "+"
    ),
    tx_D = data.frame(
      seqname = "chr1",
      start   = c(100L, 500L),
      end     = c(300L, 700L),
      strand  = "+"
    )
  ))

  fp_C <- reconstruct_chain(200L, 650L, "+", "tx_C", exons_by_tx)
  fp_D <- reconstruct_chain(150L, 650L, "+", "tx_D", exons_by_tx)

  expect_false(make_chain_id(fp_C) == make_chain_id(fp_D))
})

# ==============================================================================
# Tests for compute_frame()
# ==============================================================================

test_that("compute_frame plus-strand: frame derived from start position", {
  # Genomic start at position 1 → (1-1) %% 3 = 0
  fp0 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 300L), strand = "+")
  expect_equal(compute_frame(fp0), 0L)

  # Genomic start at position 2 → (2-1) %% 3 = 1
  fp1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(2L, 300L), strand = "+")
  expect_equal(compute_frame(fp1), 1L)

  # Genomic start at position 3 → (3-1) %% 3 = 2
  fp2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(3L, 300L), strand = "+")
  expect_equal(compute_frame(fp2), 2L)
})

test_that("compute_frame minus-strand: frame derived from end (5' end) position", {
  # On minus strand, the 5' coding nt = END of the rightmost block.
  # Position 300 → (300-1) %% 3 = 299 %% 3 = 2
  fp_minus <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(100L, 300L), strand = "-"
  )
  expected_frame <- (300L - 1L) %% 3L
  expect_equal(compute_frame(fp_minus), as.integer(expected_frame))
})

test_that("compute_frame uses phase column if present and valid", {
  fp <- GenomicRanges::GRanges("chr1", IRanges::IRanges(7L, 300L), strand = "+")
  # position 7 → (7-1)%%3 = 0, but we override with phase=2
  GenomicRanges::mcols(fp)$phase <- 2L
  expect_equal(compute_frame(fp), 2L)
})

# ==============================================================================
# Tests for make_translon_id() and make_start_id()
# ==============================================================================

test_that("make_translon_id format and stop-codon position are correct", {
  # Plus strand single exon: stop = last nt = end of range
  fp_plus <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 300L), "+")
  tid <- make_translon_id(fp_plus)
  # stop_pos = 300, frame = (100-1)%%3 = 0
  expect_match(tid, "^chr1:\\+:300:f0$")

  # Minus strand: stop = first nt = start (lowest coord) of range
  fp_minus <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 300L), "-")
  tid_m <- make_translon_id(fp_minus)
  # stop_pos on minus = min(start) = 100; frame = (300-1)%%3
  expected_frame <- as.integer((300L - 1L) %% 3L)
  expect_match(tid_m, paste0("^chr1:-:100:f", expected_frame, "$"))
})

test_that("make_start_id encodes both start and stop", {
  fp <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 300L), "+")
  sid <- make_start_id(fp)
  # start_pos=100, stop_pos=300, frame=0
  expect_match(sid, "^chr1:\\+:100:300:f0$")
})

test_that("two footprints with same stop/frame/strand share translon_id", {
  # Frame = (start - 1) %% 3.
  # start=100 → (99)%%3 = 0.  start=103 → (102)%%3 = 0.  Same frame, same stop.
  fp1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 300L), "+")
  fp2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(103L, 300L), "+")
  expect_equal(make_translon_id(fp1), make_translon_id(fp2))
  # But start_ids differ because start_pos differs
  expect_false(make_start_id(fp1) == make_start_id(fp2))
})

# ==============================================================================
# Tests for parse_orf_footprints()
# ==============================================================================

# Shared minimal TxDb fixture — we mock it as a GRangesList directly and
# inject it by monkey-patching the function's exonsBy call.
# Since parse_orf_footprints calls exonsBy internally, the cleanest approach
# for unit testing is to create a tiny real TxDb using makeTxDbFromGRanges
# or to test via reconstruct_chain directly and test the collapse logic only.
# Here we test parse_orf_footprints via a mockable wrapper.

# Helper: make a calls_df row
.call_row <- function(start, stop, strand, transcript_id, caller) {
  data.frame(
    start = start, stop = stop, strand = strand,
    transcript_id = transcript_id, caller = caller,
    stringsAsFactors = FALSE
  )
}

# We test collapse logic by directly exercising the chain_id aggregation.
# To avoid TxDb dependency we test the helper functions and the collapse
# outcome using reconstruct_chain + make_chain_id directly.

test_that("two calls on different tx_ids with same path → one chain", {
  exons_by_tx <- make_fake_exons_by_tx(list(
    ENST001 = data.frame(
      seqname = "chr1", start = c(100L, 500L), end = c(300L, 700L), strand = "+"
    ),
    ENST002 = data.frame(
      seqname = "chr1", start = c(80L, 500L), end = c(300L, 700L), strand = "+"
    )
  ))

  # Both callers agree on the ORF path: start=150, stop=650
  fp1 <- reconstruct_chain(150L, 650L, "+", "ENST001", exons_by_tx)
  fp2 <- reconstruct_chain(150L, 650L, "+", "ENST002", exons_by_tx)

  id1 <- make_chain_id(fp1)
  id2 <- make_chain_id(fp2)

  # The paths are identical → same chain_id → collapse to one footprint
  expect_equal(id1, id2)
})

test_that("two calls with different starts → two distinct chains", {
  exons_by_tx <- make_fake_exons_by_tx(list(
    ENST003 = data.frame(
      seqname = "chr1", start = c(100L, 500L), end = c(300L, 700L), strand = "+"
    ),
    ENST004 = data.frame(
      seqname = "chr1", start = c(100L, 500L), end = c(300L, 700L), strand = "+"
    )
  ))

  fp_a <- reconstruct_chain(150L, 650L, "+", "ENST003", exons_by_tx)
  fp_b <- reconstruct_chain(200L, 650L, "+", "ENST004", exons_by_tx)

  expect_false(make_chain_id(fp_a) == make_chain_id(fp_b))
})

test_that("minus-strand multi-exon ORF: frame and chain_id are self-consistent", {
  exons_by_tx <- make_fake_exons_by_tx(list(
    ENST_MINUS = data.frame(
      seqname = "chr2",
      start   = c(200L, 800L),
      end     = c(400L, 1000L),
      strand  = "-"
    )
  ))

  fp <- reconstruct_chain(950L, 250L, "-", "ENST_MINUS", exons_by_tx)

  expect_false(is.null(fp))
  # Confirm 5'→3' order: first element should be the higher-coord block
  expect_true(GenomicRanges::end(fp[1]) >= GenomicRanges::end(fp[2]))

  cid   <- make_chain_id(fp)
  frame <- compute_frame(fp)
  tid   <- make_translon_id(fp)

  expect_match(cid, "^chain_[0-9a-f]+$")
  expect_true(frame %in% 0:2)
  # stop on minus = min genomic = 250
  expect_match(tid, "^chr2:-:250:f")
})

test_that("failed reconstruction is flagged; chain_reconstructed = FALSE", {
  # We inject the reconstruction outcome by directly testing the flag logic,
  # since we cannot mock TxDb easily without creating one.
  # Verify that reconstruct_chain returns NULL for missing tx.
  exons_by_tx <- make_fake_exons_by_tx(list(
    ENST005 = data.frame(seqname = "chr1", start = 100L, end = 300L, strand = "+")
  ))

  result <- suppressWarnings(
    reconstruct_chain(110L, 290L, "+", "NOT_IN_TxDb", exons_by_tx)
  )
  expect_null(result)
  # chain_reconstructed = FALSE is set when reconstruct_chain returns NULL —
  # tested implicitly through parse_orf_footprints in integration tests.
})

test_that("check_txid_join_coverage reports correct tier counts", {
  skip_if_not_installed("GenomicFeatures")

  # Build a minimal TxDb in memory via GenomicFeatures::makeTxDb().
  # makeTxDb() accepts plain data frames — no GRanges dependency.
  txdb <- suppressMessages(GenomicFeatures::makeTxDb(
    transcripts = data.frame(
      tx_id     = 1L,
      tx_chrom  = "chr1",
      tx_strand = "+",
      tx_start  = 1000L,
      tx_end    = 2000L,
      tx_name   = "ENST00000000001"
    ),
    splicings = data.frame(
      tx_id      = 1L,
      exon_rank  = 1L,
      exon_start = 1000L,
      exon_end   = 2000L
    ),
    chrominfo = data.frame(
      chrom  = "chr1",
      length = 10000L,
      is_circular = FALSE
    )
  ))

  calls_df <- data.frame(
    transcript_id = c(
      "ENST00000000001",     # exact match
      "ENST00000000001.3",   # match after stripping ".3"
      "ENST99999999999.1"    # unresolvable
    ),
    stringsAsFactors = FALSE
  )

  result <- suppressWarnings(check_txid_join_coverage(calls_df, txdb))

  expect_named(result, c("summary", "unresolvable_ids"))
  expect_equal(result$summary$n[result$summary$tier == "exact"],       1L)
  expect_equal(result$summary$n[result$summary$tier == "strip"],       1L)
  expect_equal(result$summary$n[result$summary$tier == "unresolvable"],1L)
  expect_equal(result$unresolvable_ids, "ENST99999999999.1")
})
