make_type_reference <- function() {
  transcripts <- GenomicRanges::GRangesList(
    tx_coding = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100, 499), strand = "+"
    ),
    tx_noncoding = GenomicRanges::GRanges(
      "chr2", IRanges::IRanges(100, 499), strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_coding = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 399), strand = "+"
    )
  )
  list(transcripts = transcripts, cds = cds)
}

test_that("all supported ORF geometries are classified", {
  ref <- make_type_reference()
  orfs <- GenomicRanges::GRanges(
    seqnames = c(rep("chr1", 11), "chr2"),
    ranges = IRanges::IRanges(
      start = c(200, 197, 203, 200, 200, 197, 150, 180, 420, 350, 230, 200),
      end = c(399, 399, 399, 402, 396, 402, 180, 230, 450, 420, 350, 250)
    ),
    strand = "+"
  )
  names(orfs) <- c(
    "annotated", "n_ext", "n_trunc", "c_ext", "c_trunc", "nc_ext",
    "uorf", "uoorf", "dorf", "doorf", "intorf", "varrna"
  )

  observed <- reannotate_orf_type(
    orfs, ref$transcripts, ref$cds
  )$orfs

  expect_equal(
    as.character(observed$reference_orf_type),
    c(
      "annotated CDS",
      "N-terminal extension", "N-terminal truncation",
      "C-terminal extension", "C-terminal truncation",
      "NC-terminal extension",
      "uORF", "uoORF", "dORF", "doORF", "intORF", "varRNA-ORF"
    )
  )
  expect_identical(levels(observed$reference_orf_type), orf_type_levels())
})

test_that("reference CDS interpretation wins across transcript isoforms", {
  tx_range <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(100, 499), strand = "+"
  )
  transcripts <- GenomicRanges::GRangesList(
    tx_coding = tx_range,
    tx_noncoding = tx_range
  )
  cds <- GenomicRanges::GRangesList(
    tx_coding = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 399), strand = "+"
    )
  )
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(200, 399), strand = "+"
  )
  names(orf) <- "orf1"

  result <- reannotate_orf_type(orf, transcripts, cds)

  expect_equal(as.character(result$orfs$reference_orf_type), "annotated CDS")
  expect_equal(result$orfs$matched_transcript_id, "tx_coding")
  expect_true(result$orfs$reference_type_ambiguous)
  expect_equal(result$orfs$n_compatible_transcripts, 2L)
  expect_setequal(
    result$pairs$reference_orf_type,
    c("annotated CDS", "varRNA-ORF")
  )
})

test_that("minus-strand coordinates are interpreted in transcript direction", {
  transcripts <- GenomicRanges::GRangesList(
    tx_minus = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100, 499), strand = "-"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_minus = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 399), strand = "-"
    )
  )
  orfs <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(200, 200), end = c(399, 402)),
    strand = "-"
  )
  names(orfs) <- c("exact", "n_extension")

  observed <- reannotate_orf_type(orfs, transcripts, cds)$orfs

  expect_equal(
    as.character(observed$reference_orf_type),
    c("annotated CDS", "N-terminal extension")
  )
})

test_that("spliced transcript coordinates ignore intron widths", {
  transcripts <- GenomicRanges::GRangesList(
    tx_spliced = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(100, 300), end = c(199, 399)),
      strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_spliced = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(199, 350)),
      strand = "+"
    )
  )
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(150, 350), strand = "+"
  )
  names(orf) <- "spliced_exact"

  result <- reannotate_orf_type(orf, transcripts, cds)

  expect_equal(as.character(result$orfs$reference_orf_type), "annotated CDS")
  expect_equal(result$pairs$orf_5p_tx, result$pairs$cds_5p_tx)
  expect_equal(result$pairs$orf_3p_tx, result$pairs$cds_3p_tx)
})

test_that("an ORF without a compatible transcript receives the fallback", {
  ref <- make_type_reference()
  orf <- GenomicRanges::GRanges(
    "chr3", IRanges::IRanges(100, 160), strand = "+"
  )
  names(orf) <- "intergenic"

  result <- reannotate_orf_type(orf, ref$transcripts, ref$cds)

  expect_equal(as.character(result$orfs$reference_orf_type), "varRNA-ORF")
  expect_equal(result$orfs$n_compatible_transcripts, 0L)
  expect_match(result$orfs$reference_annotation_status, "fallback")
  expect_equal(nrow(result$pairs), 0L)
})

test_that("ORF IDs can be supplied in a metadata column", {
  ref <- make_type_reference()
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(200, 399), strand = "+",
    ORF_id = "stable_id"
  )

  result <- reannotate_orf_type(
    orf, ref$transcripts, ref$cds, id_col = "ORF_id"
  )

  expect_equal(result$orfs$orf_id, "stable_id")
})

test_that("input and priority validation fail clearly", {
  ref <- make_type_reference()
  unstranded <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(200, 399), strand = "*"
  )
  expect_error(
    reannotate_orf_type(unstranded, ref$transcripts, ref$cds),
    "explicit"
  )

  stranded <- unstranded
  GenomicRanges::strand(stranded) <- "+"
  expect_error(
    reannotate_orf_type(
      stranded, ref$transcripts, ref$cds,
      type_priority = c("annotated CDS", "varRNA-ORF")
    ),
    "every supported ORF type"
  )
})
