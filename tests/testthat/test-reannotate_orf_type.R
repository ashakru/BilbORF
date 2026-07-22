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
    orfs, ref$transcripts, ref$cds,
    stop_codon_convention = "included"
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

test_that("stop-exclusive ORFs are normalised before classification", {
  transcripts <- GenomicRanges::GRangesList(
    tx_spliced = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(100, 300), end = c(199, 499)),
      strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_spliced = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(199, 399)),
      strand = "+"
    )
  )
  orfs <- GenomicRanges::GRangesList(
    n_truncation_stop_excluded = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(153, 300), end = c(199, 396)),
      strand = "+"
    ),
    exact_stop_excluded = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(199, 396)),
      strand = "+"
    ),
    exact_stop_included = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(199, 399)),
      strand = "+"
    )
  )

  automatic <- reannotate_orf_type(orfs, transcripts, cds)

  expect_equal(
    as.character(automatic$orfs$reference_orf_type),
    c("N-terminal truncation", "annotated CDS", "annotated CDS")
  )
  expect_equal(
    automatic$orfs$reference_stop_boundary_adjusted,
    c(TRUE, TRUE, FALSE)
  )
  expect_equal(
    automatic$orfs$reference_stop_boundary_adjustment_nt,
    c(3L, 3L, 0L)
  )
  expect_equal(
    automatic$pairs$orf_3p_tx_classification,
    automatic$pairs$orf_3p_tx + c(3L, 3L, 0L)
  )
  expect_match(
    automatic$orfs$reference_annotation_status[1],
    "3 nt stop-codon adjustment"
  )

  strict <- reannotate_orf_type(
    orfs, transcripts, cds,
    stop_codon_convention = "included"
  )
  expect_equal(
    as.character(strict$orfs$reference_orf_type),
    c("intORF", "C-terminal truncation", "annotated CDS")
  )
  expect_false(any(strict$orfs$reference_stop_boundary_adjusted))
})

test_that("stop-boundary normalisation is strand aware", {
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
  orfs <- GenomicRanges::GRangesList(
    minus_n_truncation = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(203, 396), strand = "-"
    )
  )

  result <- reannotate_orf_type(orfs, transcripts, cds)

  expect_equal(
    as.character(result$orfs$reference_orf_type),
    "N-terminal truncation"
  )
  expect_true(result$orfs$reference_stop_boundary_adjusted)
  expect_equal(result$orfs$reference_stop_boundary_adjustment_nt, 3L)
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
  expect_false(result$orfs$splice_chain_supplied)
  expect_false(result$pairs$splice_chain_checked)
  expect_true(result$pairs$splice_chain_compatible)
})

test_that("exon-resolved ORFs require a compatible splice chain", {
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
  orfs <- GenomicRanges::GRangesList(
    `ORF_100000:c1` = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(199, 350)),
      strand = "+",
      exon_rank = 1:2
    ),
    `ORF_100000:c2` = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(189, 350)),
      strand = "+",
      exon_rank = 1:2
    )
  )

  result <- reannotate_orf_type(orfs, transcripts, cds)

  expect_equal(
    result$orfs$orf_id,
    c("ORF_100000:c1", "ORF_100000:c2")
  )
  expect_equal(result$orfs$n_orf_exons, c(2L, 2L))
  expect_true(all(result$orfs$splice_chain_supplied))
  expect_equal(
    as.character(result$orfs$reference_orf_type),
    c("annotated CDS", NA_character_)
  )
  expect_equal(result$orfs$n_compatible_transcripts, c(1L, 0L))

  compatible <- result$pairs[result$pairs$orf_id == "ORF_100000:c1", ]
  incompatible <- result$pairs[result$pairs$orf_id == "ORF_100000:c2", ]
  expect_true(compatible$splice_chain_checked)
  expect_true(compatible$splice_chain_compatible)
  expect_true(incompatible$splice_chain_checked)
  expect_false(incompatible$splice_chain_compatible)
  expect_match(incompatible$annotation_status, "splice chain")
  expect_match(
    result$orfs$reference_annotation_status[2],
    "no transcript has a compatible ORF splice chain"
  )
})

test_that("exon-resolved minus-strand ORFs retain transcript orientation", {
  transcripts <- GenomicRanges::GRangesList(
    tx_minus = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(100, 300), end = c(199, 399)),
      strand = "-"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_minus = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(199, 350)),
      strand = "-"
    )
  )
  orfs <- GenomicRanges::GRangesList(
    minus_exact = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(199, 350)),
      strand = "-",
      exon_rank = 2:1
    )
  )

  result <- reannotate_orf_type(orfs, transcripts, cds)

  expect_equal(
    as.character(result$orfs$reference_orf_type),
    "annotated CDS"
  )
  expect_true(result$pairs$splice_chain_checked)
  expect_true(result$pairs$splice_chain_compatible)
  expect_equal(result$pairs$orf_5p_tx, result$pairs$cds_5p_tx)
  expect_equal(result$pairs$orf_3p_tx, result$pairs$cds_3p_tx)
})

test_that("short terminal overhangs infer canonical extensions", {
  transcripts <- GenomicRanges::GRangesList(
    tx_short = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100, 499), strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_short = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 399), strand = "+"
    )
  )
  orfs <- GenomicRanges::GRangesList(
    n_extension = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(95, 399), strand = "+"
    ),
    c_extension = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 504), strand = "+"
    ),
    nc_extension = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(95, 504), strand = "+"
    )
  )

  result <- reannotate_orf_type(orfs, transcripts, cds)

  expect_equal(
    as.character(result$orfs$reference_orf_type),
    c(
      "N-terminal extension", "C-terminal extension",
      "NC-terminal extension"
    )
  )
  expect_true(all(result$orfs$reference_boundary_extrapolated))
  expect_equal(result$orfs$reference_5p_overhang_nt, c(5L, 0L, 5L))
  expect_equal(result$orfs$reference_3p_overhang_nt, c(0L, 5L, 5L))
  expect_true(all(result$pairs$terminal_extrapolation_used))
  expect_true(all(grepl(
    "classified using terminal extrapolation",
    result$orfs$reference_annotation_status
  )))
})

test_that("terminal extrapolation preserves a compatible splice chain", {
  transcripts <- GenomicRanges::GRangesList(
    tx_spliced = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(100, 300), end = c(199, 499)),
      strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_spliced = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(150, 300), end = c(199, 399)),
      strand = "+"
    )
  )
  orfs <- GenomicRanges::GRangesList(
    spliced_extension = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = c(99, 300), end = c(199, 399)),
      strand = "+",
      exon_rank = 1:2
    )
  )

  result <- reannotate_orf_type(orfs, transcripts, cds)

  expect_equal(
    as.character(result$orfs$reference_orf_type),
    "N-terminal extension"
  )
  expect_true(result$orfs$reference_boundary_extrapolated)
  expect_true(result$orfs$reference_splice_chain_compatible)
  expect_equal(result$orfs$reference_5p_overhang_nt, 1L)
})

test_that("terminal extrapolation respects minus-strand orientation", {
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
  orfs <- GenomicRanges::GRangesList(
    minus_n_extension = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 504), strand = "-"
    )
  )

  result <- reannotate_orf_type(orfs, transcripts, cds)

  expect_equal(
    as.character(result$orfs$reference_orf_type),
    "N-terminal extension"
  )
  expect_true(result$orfs$reference_boundary_extrapolated)
  expect_equal(result$orfs$reference_5p_overhang_nt, 5L)
  expect_equal(result$orfs$reference_3p_overhang_nt, 0L)
})

test_that("a complete longer transcript takes precedence over extrapolation", {
  transcripts <- GenomicRanges::GRangesList(
    tx_short = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100, 499), strand = "+"
    ),
    tx_long = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(90, 499), strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_short = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 399), strand = "+"
    ),
    tx_long = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 399), strand = "+"
    )
  )
  orfs <- GenomicRanges::GRangesList(
    candidate = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(95, 399), strand = "+"
    )
  )

  result <- reannotate_orf_type(orfs, transcripts, cds)

  expect_equal(result$orfs$matched_transcript_id, "tx_long")
  expect_equal(
    as.character(result$orfs$reference_orf_type),
    "N-terminal extension"
  )
  expect_false(result$orfs$reference_boundary_extrapolated)
  expect_equal(result$orfs$n_compatible_transcripts, 1L)
  expect_true(
    result$pairs$terminal_extrapolation_used[
      result$pairs$transcript_id == "tx_short"
    ]
  )
})

test_that("invalid terminal extrapolations remain unclassified", {
  transcripts <- GenomicRanges::GRangesList(
    tx_short = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(100, 499), strand = "+"
    )
  )
  cds <- GenomicRanges::GRangesList(
    tx_short = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(200, 399), strand = "+"
    )
  )
  orfs <- GenomicRanges::GRangesList(
    over_limit = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(65, 399), strand = "+"
    ),
    out_of_frame = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(96, 399), strand = "+"
    )
  )

  result <- reannotate_orf_type(orfs, transcripts, cds)

  expect_true(all(is.na(result$orfs$reference_orf_type)))
  expect_true(all(is.na(result$orfs$reference_orf_class)))
  expect_false(any(result$pairs$terminal_extrapolation_used))
  expect_match(result$pairs$annotation_status[1], "exceeds")
  expect_match(result$pairs$annotation_status[2], "in-frame canonical")

  disabled <- reannotate_orf_type(
    orfs[2], transcripts, cds,
    allow_terminal_extrapolation = FALSE
  )
  expect_true(is.na(disabled$orfs$reference_orf_type))
  expect_match(disabled$pairs$annotation_status, "disabled")
})

test_that("an ORF without a compatible transcript remains unclassified", {
  ref <- make_type_reference()
  orf <- GenomicRanges::GRanges(
    "chr3", IRanges::IRanges(100, 160), strand = "+"
  )
  names(orf) <- "intergenic"

  result <- reannotate_orf_type(orf, ref$transcripts, ref$cds)

  expect_true(is.na(result$orfs$reference_orf_type))
  expect_true(is.na(result$orfs$reference_orf_class))
  expect_equal(result$orfs$n_compatible_transcripts, 0L)
  expect_match(result$orfs$reference_annotation_status, "unclassified")
  expect_equal(nrow(result$pairs), 0L)

  legacy <- reannotate_orf_type(
    orf, ref$transcripts, ref$cds, unmatched_type = "varRNA-ORF"
  )
  expect_equal(as.character(legacy$orfs$reference_orf_type), "varRNA-ORF")
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

test_that("unused ORF seqlevels do not break transcript mapping", {
  ref <- make_type_reference()
  orf <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(200, 399), strand = "+"
  )
  GenomeInfoDb::seqlevels(orf) <- c("chr1", "unused_contig")
  names(orf) <- "extra_seqlevel"

  result <- reannotate_orf_type(orf, ref$transcripts, ref$cds)

  expect_equal(
    as.character(result$orfs$reference_orf_type),
    "annotated CDS"
  )
  expect_equal(result$orfs$matched_transcript_id, "tx_coding")
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

  mixed_seqnames <- GenomicRanges::GRangesList(
    bad_orf = GenomicRanges::GRanges(
      c("chr1", "chr2"),
      IRanges::IRanges(c(200, 300), c(250, 350)),
      strand = "+"
    )
  )
  expect_error(
    reannotate_orf_type(mixed_seqnames, ref$transcripts, ref$cds),
    "one seqname and one strand"
  )
  expect_error(
    reannotate_orf_type(
      stranded, ref$transcripts, ref$cds,
      max_terminal_overhang = -1L
    ),
    "non-negative integer"
  )
  expect_error(
    reannotate_orf_type(
      stranded, ref$transcripts, ref$cds,
      unmatched_type = "unsupported"
    ),
    "supported ORF type"
  )
  expect_error(
    reannotate_orf_type(
      stranded, ref$transcripts, ref$cds,
      stop_codon_convention = "unknown"
    ),
    "arg"
  )
})
