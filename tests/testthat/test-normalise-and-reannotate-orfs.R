test_that("stop normalisation reuses isoform geometry across a splice junction", {
  genome_chars <- rep("A", 300L)
  genome_chars[c(101L, 102L, 103L, 104L, 200L, 201L)] <-
    strsplit("ATGTGA", "", fixed = TRUE)[[1L]]
  genome <- Biostrings::DNAStringSet(paste0(genome_chars, collapse = ""))
  names(genome) <- "chr1"

  transcripts <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(100L, 200L), c(104L, 209L)),
      strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(101L, 200L), c(104L, 201L)),
      strand = "+"
    )
  )
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(101L, 103L), strand = "+"
  )
  names(orf) <- "orf1"

  observed <- normalise_and_reannotate_orfs(
    orfs = orf,
    transcripts = transcripts,
    cds_by_tx = cds,
    genome = genome,
    stop_codon_convention = "excluded"
  )

  expect_equal(observed$orfs$canonical_stop_triplet, "TGA")
  expect_true(observed$orfs$canonical_stop_valid)
  expect_true(observed$orfs$canonical_complete_codons)
  expect_equal(
    as.character(observed$orfs$reference_orf_type), "annotated CDS"
  )
  expect_equal(length(observed$canonical_chains), 1L)
  expect_equal(unname(lengths(observed$canonical_chains)), 2L)
  expect_equal(
    sum(GenomicRanges::width(observed$canonical_chains[[1L]])), 6L
  )
  expect_equal(
    unname(as.character(GenomicFeatures::extractTranscriptSeqs(
      genome, observed$canonical_chains
    ))),
    "ATGTGA"
  )
})

test_that("invalid following triplets do not produce canonical chains", {
  genome <- Biostrings::DNAStringSet(paste(rep("A", 300L), collapse = ""))
  names(genome) <- "chr1"
  transcripts <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100L, 299L), strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100L, 108L), strand = "+"
    )
  )
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(100L, 105L), strand = "+"
  )
  names(orf) <- "orf1"

  observed <- normalise_and_reannotate_orfs(
    orf, transcripts, cds, genome,
    stop_codon_convention = "excluded"
  )

  expect_false(observed$orfs$canonical_stop_valid)
  expect_length(observed$canonical_chains, 0L)
  expect_match(
    observed$orfs$reference_annotation_status,
    "valid canonical stop"
  )
})

test_that("canonical chain construction is strand aware", {
  genome_chars <- rep("A", 300L)
  genome_chars[101L:109L] <- strsplit("TCATTTCAT", "", fixed = TRUE)[[1L]]
  genome <- Biostrings::DNAStringSet(paste0(genome_chars, collapse = ""))
  names(genome) <- "chr1"

  transcripts <- GenomicRanges::GRangesList(
    tx_minus = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100L, 109L), strand = "-"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_minus = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(101L, 109L), strand = "-"
    )
  )
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(104L, 109L), strand = "-"
  )
  names(orf) <- "orf_minus"

  observed <- normalise_and_reannotate_orfs(
    orf, transcripts, cds, genome,
    stop_codon_convention = "excluded"
  )

  expect_equal(observed$orfs$canonical_stop_triplet, "TGA")
  expect_equal(
    unname(as.character(GenomicFeatures::extractTranscriptSeqs(
      genome, observed$canonical_chains
    ))),
    "ATGAAATGA"
  )
  expect_equal(GenomicRanges::start(observed$canonical_chains[[1L]]), 101L)
  expect_equal(GenomicRanges::end(observed$canonical_chains[[1L]]), 109L)
})
