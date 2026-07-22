make_stop_coordinate_fixture <- function() {
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
  list(
    genome = genome,
    annotations = list(transcripts = transcripts),
    cds = GenomicRanges::GRangesList(
      tx1 = GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(101L, 200L), c(104L, 201L)),
        strand = "+"
      )
    ),
    transcript_meta = data.frame(
      transcript_id = "tx1", gene_id = "g1", gene_name = "G1"
    )
  )
}

test_that("annotate_orf_isoforms returns stop-inclusive ranges", {
  fixture <- make_stop_coordinate_fixture()
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(101L, 103L), strand = "+"
  )
  names(orf) <- "orf1"

  observed <- annotate_orf_isoforms(
    annotations = fixture$annotations,
    orfs = orf,
    BSgenome = fixture$genome,
    transcript_meta = fixture$transcript_meta,
    cds_gr = fixture$cds,
    stop_codon_convention = "excluded"
  )

  expect_equal(length(observed$ranges), 1L)
  expect_equal(unname(lengths(observed$ranges)), 2L)
  expect_equal(
    unname(as.character(GenomicFeatures::extractTranscriptSeqs(
      fixture$genome, observed$ranges
    ))),
    "ATGTGA"
  )
  expect_equal(observed$table$input_stop_convention, "excluded")
  expect_equal(observed$table$stop_adjustment_nt, 3L)
  expect_true(observed$table$stop_normalization_valid)
  expect_equal(observed$table$canonical_stop_codon, "TGA")
  expect_equal(observed$table$start_stop_exclusive, 101L)
  expect_equal(observed$table$end_stop_exclusive, 103L)
  expect_equal(observed$table$orf_3p_stop_exclusive, 103L)
  expect_equal(observed$table$start_stop_inclusive, 101L)
  expect_equal(observed$table$end_stop_inclusive, 201L)
  expect_equal(observed$table$orf_3p_stop_inclusive, 201L)
  expect_equal(observed$table$start, observed$table$start_stop_inclusive)
  expect_equal(observed$table$end, observed$table$end_stop_inclusive)
  expect_equal(observed$table$canonical_5p_tx, 2L)
  expect_equal(observed$table$canonical_3p_tx, 7L)
  expect_equal(as.character(observed$table$reference_orf_type), "annotated CDS")
  expect_equal(as.character(observed$table$reference_orf_class), "canonical")
})

test_that("included and auto conventions report both coordinate systems", {
  fixture <- make_stop_coordinate_fixture()
  orf <- GenomicRanges::GRangesList(
    orf1 = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(101L, 200L), c(104L, 201L)),
      strand = "+"
    )
  )

  included <- annotate_orf_isoforms(
    fixture$annotations, orf, fixture$genome, fixture$transcript_meta,
    stop_codon_convention = "included"
  )
  automatic <- annotate_orf_isoforms(
    fixture$annotations, orf, fixture$genome, fixture$transcript_meta,
    stop_codon_convention = "auto"
  )

  expect_equal(included$table$input_stop_convention, "included")
  expect_equal(included$table$stop_adjustment_nt, 0L)
  expect_equal(included$table$end_stop_exclusive, 103L)
  expect_equal(included$table$end_stop_inclusive, 201L)
  expect_equal(automatic$table$input_stop_convention, "included")
  expect_true(automatic$table$stop_normalization_valid)
})

test_that("auto resolves a mixture record by record", {
  fixture <- make_stop_coordinate_fixture()
  mixed <- GenomicRanges::GRangesList(
    excluded = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(101L, 103L), strand = "+"
    ),
    included = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(101L, 200L), c(104L, 201L)),
      strand = "+"
    )
  )

  observed <- annotate_orf_isoforms(
    fixture$annotations, mixed, fixture$genome, fixture$transcript_meta,
    stop_codon_convention = "auto"
  )

  convention <- stats::setNames(
    observed$table$input_stop_convention, observed$table$ORF_id
  )
  expect_equal(unname(convention[c("excluded", "included")]),
               c("excluded", "included"))
  expect_true(all(observed$table$stop_normalization_valid))
  expect_equal(
    unname(as.character(GenomicFeatures::extractTranscriptSeqs(
      fixture$genome, observed$ranges
    ))),
    c("ATGTGA", "ATGTGA")
  )
})

test_that("all compatible transcript isoforms are retained", {
  fixture <- make_stop_coordinate_fixture()
  fixture$annotations$transcripts <- GenomicRanges::GRangesList(
    tx1 = fixture$annotations$transcripts[["tx1"]],
    tx2 = fixture$annotations$transcripts[["tx1"]]
  )
  fixture$cds <- GenomicRanges::GRangesList(
    tx1 = fixture$cds[["tx1"]],
    tx2 = fixture$cds[["tx1"]]
  )
  fixture$transcript_meta <- data.frame(
    transcript_id = c("tx1", "tx2"),
    gene_id = c("g1", "g1"),
    gene_name = c("G1", "G1")
  )
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(101L, 103L), strand = "+"
  )
  names(orf) <- "orf1"

  observed <- annotate_orf_isoforms(
    fixture$annotations, orf, fixture$genome, fixture$transcript_meta,
    cds_gr = fixture$cds,
    stop_codon_convention = "auto"
  )

  expect_setequal(observed$table$transcript_id, c("tx1", "tx2"))
  expect_equal(length(observed$ranges), 2L)
  expect_true(all(observed$table$stop_normalization_valid))
  expect_true(all(as.character(observed$table$reference_orf_type) ==
                    "annotated CDS"))
})

test_that("stop normalization is transcript-oriented on the minus strand", {
  genome_chars <- rep("A", 300L)
  # Reverse complement of transcript-oriented ATGAAATGA.
  genome_chars[101L:109L] <-
    strsplit("TCATTTCAT", "", fixed = TRUE)[[1L]]
  genome <- Biostrings::DNAStringSet(paste0(genome_chars, collapse = ""))
  names(genome) <- "chr1"
  transcripts <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100L, 109L), strand = "-"
    )
  )
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(104L, 109L), strand = "-"
  )
  names(orf) <- "orf1"

  observed <- annotate_orf_isoforms(
    annotations = list(transcripts = transcripts),
    orfs = orf,
    BSgenome = genome,
    transcript_meta = data.frame(
      transcript_id = "tx1", gene_id = "g1", gene_name = "G1"
    ),
    stop_codon_convention = "excluded"
  )

  expect_equal(
    unname(as.character(GenomicFeatures::extractTranscriptSeqs(
      genome, observed$ranges
    ))),
    "ATGAAATGA"
  )
  expect_equal(observed$table$start_stop_exclusive, 104L)
  expect_equal(observed$table$end_stop_exclusive, 109L)
  expect_equal(observed$table$orf_3p_stop_exclusive, 104L)
  expect_equal(observed$table$start_stop_inclusive, 101L)
  expect_equal(observed$table$end_stop_inclusive, 109L)
  expect_equal(observed$table$orf_3p_stop_inclusive, 101L)
  expect_true(observed$table$stop_normalization_valid)
})
