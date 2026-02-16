# ==============================================================================
# Tests for the spec-based ORF parser system
# ==============================================================================

library(testthat)
library(GenomicRanges)

# -- Mock data helpers ---------------------------------------------------------

create_ribotish_test_file <- function(dir) {
  file <- file.path(dir, "ribotish.txt")
  df <- data.frame(
    GenomePos   = c("chr1:100-300:+", "chr1:500-800:-", "chr2:200-400:+"),
    Tid         = c("ENST001", "ENST002", "ENST003"),
    Symbol      = c("GENE1", "GENE2", "GENE3"),
    TisType     = c("5'UTR", "CDS", "3'UTR"),
    StartCodon  = c("ATG", "ATG", "CTG"),
    AALen       = c(66, 100, 66),
    RiboPvalue  = c(0.001, 0.01, 0.05),
    stringsAsFactors = FALSE
  )
  write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
  file
}

create_ribocode_test_file <- function(dir) {
  file <- file.path(dir, "ribocode.txt")
  df <- data.frame(
    ORF_ID        = c("orf_1", "orf_2", "orf_3"),
    chrom         = c("chr1", "chr1", "chr2"),
    ORF_gstart    = c(100, 500, 200),
    ORF_gstop     = c(300, 800, 400),
    strand        = c("+", "-", "+"),
    ORF_type      = c("uORF", "CDS", "dORF"),
    transcript_id = c("ENST001", "ENST002", "ENST003"),
    gene_id       = c("ENSG001", "ENSG002", "ENSG003"),
    ORF_length    = c(201, 301, 201),
    pval          = c(0.001, 0.01, 0.05),
    stringsAsFactors = FALSE
  )
  write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
  file
}

create_ribotricer_test_file <- function(dir) {
  file <- file.path(dir, "ribotricer.tsv")
  df <- data.frame(
    ORF_ID      = c("ENST001__chr1__100__300__+__uORF",
                    "ENST002__chr1__500__800__-__CDS",
                    "ENST003__chr2__200__400__+__dORF"),
    phase_score = c(0.8, 0.95, 0.7),
    stringsAsFactors = FALSE
  )
  write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
  file
}

create_orfquant_test_file <- function(dir) {
  file <- file.path(dir, "orfquant.tsv")
  df <- data.frame(
    ORF_id_tr      = c("orf_1", "orf_2", "orf_3"),
    Chromosome     = c("chr1", "chr1", "chr2"),
    ORF_Start      = c(100, 500, 200),
    ORF_End        = c(300, 800, 400),
    Strand         = c("+", "-", "+"),
    ORF_category   = c("uORF", "CDS", "dORF"),
    P_sites_raw    = c(10, 50, 5),
    stringsAsFactors = FALSE
  )
  write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
  file
}

create_bed_test_file <- function(dir) {
  file <- file.path(dir, "orfs.bed")
  lines <- c(
    "chr1\t99\t300\torf_1\t100\t+",
    "chr1\t499\t800\torf_2\t200\t-",
    "chr2\t199\t400\torf_3\t150\t+"
  )
  writeLines(lines, file)
  file
}


# ==============================================================================
# 1. orf_caller_spec construction
# ==============================================================================

test_that("orf_caller_spec creates valid S3 object", {
  spec <- orf_caller_spec(
    name       = "test_caller",
    file_format = "tsv",
    column_map = list(
      chrom = "chr", start = "s", end = "e",
      strand = "str", orf_id = "id"
    )
  )
  expect_s3_class(spec, "orf_caller_spec")
  expect_equal(spec$name, "test_caller")
  expect_equal(spec$file_format, "tsv")
  expect_equal(spec$coord_system, "1-based")
})

test_that("orf_caller_spec lowercases name", {
  spec <- orf_caller_spec(
    name = "MyTool", file_format = "bed", column_map = list()
  )
  expect_equal(spec$name, "mytool")
})

test_that("orf_caller_spec validates name", {
  expect_error(orf_caller_spec(name = ""), "non-empty")
  expect_error(orf_caller_spec(name = 42), "non-empty")
})

test_that("orf_caller_spec requires column_map fields for TSV without read_fn", {
  expect_error(
    orf_caller_spec(name = "bad", file_format = "tsv",
                    column_map = list(chrom = "chr")),
    "missing required"
  )
})

test_that("orf_caller_spec skips column validation for BED/RDS", {
  expect_no_error(
    orf_caller_spec(name = "bed_tool", file_format = "bed")
  )
  expect_no_error(
    orf_caller_spec(name = "rds_tool", file_format = "rds")
  )
})

test_that("orf_caller_spec skips column validation when read_fn is provided", {
  expect_no_error(
    orf_caller_spec(name = "custom", file_format = "tsv",
                    read_fn = function(f) data.frame())
  )
})

test_that("orf_caller_spec validates callable arguments", {
  expect_error(
    orf_caller_spec(name = "x", file_format = "bed", read_fn = "not_a_fn"),
    "function"
  )
  expect_error(
    orf_caller_spec(name = "x", file_format = "bed", post_process_fn = 42),
    "function"
  )
})

test_that("print.orf_caller_spec works", {
  spec <- orf_caller_spec(
    name = "test", file_format = "tsv",
    column_map = list(chrom = "c", start = "s", end = "e",
                      strand = "str", orf_id = "id"),
    url = "https://example.com",
    read_fn = identity,
    post_process_fn = identity
  )
  out <- capture.output(print(spec))
  expect_true(any(grepl("test", out)))
  expect_true(any(grepl("Custom reader", out)))
  expect_true(any(grepl("Post-process", out)))
})


# ==============================================================================
# 2. Registry functions
# ==============================================================================

test_that("register_orf_caller validates input", {
  expect_error(register_orf_caller("not_a_spec"), "orf_caller_spec")
})

test_that("register/get round-trip works", {
  spec <- orf_caller_spec(
    name = "roundtrip_test", file_format = "bed"
  )
  register_orf_caller(spec, overwrite = TRUE)
  retrieved <- get_orf_caller("roundtrip_test")
  expect_s3_class(retrieved, "orf_caller_spec")
  expect_equal(retrieved$name, "roundtrip_test")
})

test_that("get_orf_caller is case-insensitive", {
  spec <- orf_caller_spec(name = "casetest", file_format = "bed")
  register_orf_caller(spec, overwrite = TRUE)
  expect_s3_class(get_orf_caller("CaseTest"), "orf_caller_spec")
})

test_that("get_orf_caller errors on unknown caller", {
  expect_error(get_orf_caller("nonexistent_xyzzy"), "Unknown ORF caller")
})

test_that("register_orf_caller warns on overwrite without flag", {
  spec <- orf_caller_spec(name = "overwarn", file_format = "bed")
  register_orf_caller(spec, overwrite = TRUE)
  expect_warning(register_orf_caller(spec), "Overwriting")
})

test_that("supported_orf_callers returns expected built-in callers", {
  # Make sure builtins are registered
  .register_builtin_specs()
  callers <- supported_orf_callers()
  expect_true(is.data.frame(callers))
  expect_true(all(c("source", "description", "file_format") %in% colnames(callers)))
  expected <- c("ribotish", "ribocode", "price", "ribotricer",
                "orfquant", "gencode", "ribotie")
  expect_true(all(expected %in% callers$source))
})


# ==============================================================================
# 3. Internal helpers
# ==============================================================================

test_that(".resolve_column finds first alias match", {
  cmap <- list(chrom = c("Chr", "chrom", "seqnames"))
  available <- c("id", "chrom", "start")
  expect_equal(.resolve_column("chrom", cmap, available), "chrom")
})

test_that(".resolve_column returns NA for optional missing field", {
  cmap <- list()
  expect_true(is.na(.resolve_column("gene_id", cmap, c("a", "b"))))
})

test_that(".resolve_column errors on required missing field", {
  cmap <- list(chrom = "Chr")
  expect_error(
    .resolve_column("chrom", cmap, c("x", "y"), required = TRUE),
    "Cannot find column"
  )
})

test_that(".resolve_column errors when field not in column_map and required", {
  expect_error(
    .resolve_column("missing_field", list(), c("a", "b"), required = TRUE),
    "no entry"
  )
})


# ==============================================================================
# 4. parse_orfs — input validation
# ==============================================================================

test_that("parse_orfs errors on missing file", {
  expect_error(parse_orfs("/nonexistent.bed", "gencode"), "File not found")
})

test_that("parse_orfs errors on invalid genome_style", {
  tmp <- tempfile()
  writeLines("", tmp)
  expect_error(parse_orfs(tmp, "gencode", genome_style = "NCBI"), "genome_style")
  unlink(tmp)
})

test_that("parse_orfs errors on unknown source", {
  tmp <- tempfile()
  writeLines("", tmp)
  expect_error(parse_orfs(tmp, "totally_unknown_caller"), "Unknown ORF caller")
  unlink(tmp)
})


# ==============================================================================
# 5. parse_orfs — caller-specific parsing
# ==============================================================================

test_that("parse_orfs handles Ribo-TISH correctly", {
  dir <- tempdir()
  file <- create_ribotish_test_file(dir)

  gr <- parse_orfs(file, source = "ribotish")

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
  expect_true("orf_id" %in% colnames(mcols(gr)))
  expect_true("orf_type" %in% colnames(mcols(gr)))
  expect_true("gene_name" %in% colnames(mcols(gr)))
  expect_equal(as.character(mcols(gr)$orf_type), c("5'UTR", "CDS", "3'UTR"))

  unlink(file)
})

test_that("Ribo-TISH retains spec extra_cols", {
  dir <- tempdir()
  file <- create_ribotish_test_file(dir)

  gr <- parse_orfs(file, source = "ribotish")
  expect_true("AALen" %in% colnames(mcols(gr)))
  expect_true("StartCodon" %in% colnames(mcols(gr)))

  unlink(file)
})

test_that("parse_orfs handles RiboCode correctly", {
  dir <- tempdir()
  file <- create_ribocode_test_file(dir)

  gr <- parse_orfs(file, source = "ribocode")

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
  expect_true("orf_type" %in% colnames(mcols(gr)))
  expect_true("gene_id" %in% colnames(mcols(gr)))
  expect_true("transcript_id" %in% colnames(mcols(gr)))
  expect_equal(as.character(mcols(gr)$orf_type), c("uORF", "CDS", "dORF"))

  unlink(file)
})

test_that("parse_orfs handles PRICE (BED) correctly", {
  dir <- tempdir()
  file <- create_bed_test_file(dir)

  gr <- parse_orfs(file, source = "price")

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
  expect_true("orf_id" %in% colnames(mcols(gr)))

  unlink(file)
})

test_that("parse_orfs handles Ribotricer with parsed ORF_ID", {
  dir <- tempdir()
  file <- create_ribotricer_test_file(dir)

  gr <- parse_orfs(file, source = "ribotricer")

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
  expect_true("orf_type" %in% colnames(mcols(gr)))
  # Coordinates should have been extracted from the ORF_ID
  expect_equal(start(gr)[1], 100)
  expect_equal(end(gr)[1], 300)

  unlink(file)
})

test_that("parse_orfs handles ORFquant correctly", {
  dir <- tempdir()
  file <- create_orfquant_test_file(dir)

  gr <- parse_orfs(file, source = "orfquant")

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
  expect_true("orf_type" %in% colnames(mcols(gr)))
  expect_true("P_sites_raw" %in% colnames(mcols(gr)))

  unlink(file)
})

test_that("parse_orfs handles GENCODE BED correctly", {
  dir <- tempdir()
  file <- create_bed_test_file(dir)

  gr <- parse_orfs(file, source = "gencode")

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
  # Post-process should add orf_type

  expect_true("orf_type" %in% colnames(mcols(gr)))
  expect_true(all(mcols(gr)$orf_type == "ribo-seq_orf"))

  unlink(file)
})


# ==============================================================================
# 6. parse_orfs — additional_cols & min_length
# ==============================================================================

test_that("additional_cols are retained in metadata", {
  dir <- tempdir()
  file <- create_ribocode_test_file(dir)

  gr <- parse_orfs(file, source = "ribocode",
                   additional_cols = c("ORF_length", "pval"))
  expect_true("ORF_length" %in% colnames(mcols(gr)))
  expect_true("pval" %in% colnames(mcols(gr)))

  unlink(file)
})

test_that("additional_cols warns on missing column", {
  dir <- tempdir()
  file <- create_ribocode_test_file(dir)

  expect_warning(
    parse_orfs(file, source = "ribocode",
               additional_cols = c("nonexistent_column")),
    "not found"
  )

  unlink(file)
})

test_that("min_length filters short ORFs", {
  dir <- tempdir()
  file <- create_ribocode_test_file(dir)

  gr <- parse_orfs(file, source = "ribocode", min_length = 250)
  expect_true(all(width(gr) >= 250))

  unlink(file)
})

test_that("min_length warns when all ORFs filtered", {
  dir <- tempdir()
  file <- create_ribocode_test_file(dir)

  expect_warning(
    parse_orfs(file, source = "ribocode", min_length = 100000),
    "All ORFs filtered"
  )

  unlink(file)
})


# ==============================================================================
# 7. Custom caller registration end-to-end
# ==============================================================================

test_that("custom caller can be registered and used", {
  dir <- tempdir()

  # Create a custom TSV
  custom_file <- file.path(dir, "custom_orfs.tsv")
  df <- data.frame(
    chr    = c("chr1", "chr2"),
    begin  = c(100, 200),
    stop   = c(300, 400),
    str    = c("+", "-"),
    name   = c("orf_a", "orf_b"),
    score  = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )
  write.table(df, custom_file, sep = "\t", row.names = FALSE, quote = FALSE)

  # Define and register
  spec <- orf_caller_spec(
    name        = "test_custom_e2e",
    description = "Custom test caller",
    file_format = "tsv",
    column_map  = list(
      chrom  = "chr",
      start  = "begin",
      end    = "stop",
      strand = "str",
      orf_id = "name"
    ),
    extra_cols = "score"
  )
  register_orf_caller(spec, overwrite = TRUE)

  # Parse
  gr <- parse_orfs(custom_file, source = "test_custom_e2e")

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(as.character(names(gr)), c("orf_a", "orf_b"))
  expect_true("score" %in% colnames(mcols(gr)))

  unlink(custom_file)
})

test_that("custom caller with 0-based coords adjusts starts", {
  dir <- tempdir()

  custom_file <- file.path(dir, "zerobased.tsv")
  df <- data.frame(
    chrom  = c("chr1"),
    start  = c(99),
    end    = c(300),
    strand = c("+"),
    id     = c("orf_z"),
    stringsAsFactors = FALSE
  )
  write.table(df, custom_file, sep = "\t", row.names = FALSE, quote = FALSE)

  spec <- orf_caller_spec(
    name         = "test_zerobased",
    file_format  = "tsv",
    column_map   = list(chrom = "chrom", start = "start", end = "end",
                        strand = "strand", orf_id = "id"),
    coord_system = "0-based"
  )
  register_orf_caller(spec, overwrite = TRUE)

  gr <- parse_orfs(custom_file, source = "test_zerobased")
  expect_equal(start(gr), 100)  # 99 + 1

  unlink(custom_file)
})

test_that("custom caller with read_fn works", {
  dir <- tempdir()

  strange_file <- file.path(dir, "strange.txt")
  writeLines(c("chr1|100|300|+|orf_x", "chr2|500|700|-|orf_y"), strange_file)

  spec <- orf_caller_spec(
    name        = "test_strange_reader",
    file_format = "tsv",
    column_map  = list(
      chrom = "chrom", start = "start", end = "end",
      strand = "strand", orf_id = "id"
    ),
    read_fn = function(file) {
      lines <- readLines(file)
      parts <- strsplit(lines, "\\|")
      data.frame(
        chrom  = vapply(parts, `[`, character(1), 1),
        start  = as.integer(vapply(parts, `[`, character(1), 2)),
        end    = as.integer(vapply(parts, `[`, character(1), 3)),
        strand = vapply(parts, `[`, character(1), 4),
        id     = vapply(parts, `[`, character(1), 5),
        stringsAsFactors = FALSE
      )
    }
  )
  register_orf_caller(spec, overwrite = TRUE)

  gr <- parse_orfs(strange_file, source = "test_strange_reader")
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(as.character(names(gr)), c("orf_x", "orf_y"))

  unlink(strange_file)
})

test_that("custom caller with post_process_fn works", {
  dir <- tempdir()
  file <- create_ribocode_test_file(dir)

  spec <- orf_caller_spec(
    name        = "test_postprocess",
    file_format = "tsv",
    column_map  = list(
      chrom = "chrom", start = "ORF_gstart", end = "ORF_gstop",
      strand = "strand", orf_id = "ORF_ID"
    ),
    post_process_fn = function(gr, raw_df) {
      GenomicRanges::mcols(gr)$custom_flag <- TRUE
      gr
    }
  )
  register_orf_caller(spec, overwrite = TRUE)

  gr <- parse_orfs(file, source = "test_postprocess")
  expect_true("custom_flag" %in% colnames(mcols(gr)))
  expect_true(all(mcols(gr)$custom_flag))

  unlink(file)
})

test_that("overriding a built-in spec works", {
  original <- get_orf_caller("gencode")

  custom_gencode <- orf_caller_spec(
    name        = "gencode",
    description = "Modified GENCODE spec",
    file_format = "bed",
    column_map  = list(orf_id = "name"),
    post_process_fn = function(gr, raw_df) {
      GenomicRanges::mcols(gr)$orf_type <- "custom_type"
      gr
    }
  )
  register_orf_caller(custom_gencode, overwrite = TRUE)

  dir <- tempdir()
  file <- create_bed_test_file(dir)
  gr <- parse_orfs(file, source = "gencode")
  expect_true(all(mcols(gr)$orf_type == "custom_type"))

  # Restore original
  register_orf_caller(original, overwrite = TRUE)

  unlink(file)
})


# ==============================================================================
# 8. Column alias resolution
# ==============================================================================

test_that("column aliases resolve correctly for ORFquant", {
  dir <- tempdir()
  file <- create_orfquant_test_file(dir)

  gr <- parse_orfs(file, source = "orfquant")

  # ORFquant spec has aliases like c("Chromosome", "chromosome", ...)
  # Our test file uses "Chromosome" -- it should resolve properly
  expect_equal(as.character(seqnames(gr))[1], "chr1")

  unlink(file)
})


# ==============================================================================
# 9. ORFs are named and unique
# ==============================================================================

test_that("duplicate ORF names are made unique", {
  dir <- tempdir()
  file <- file.path(dir, "dup_ids.tsv")
  df <- data.frame(
    chrom = c("chr1", "chr1"), start = c(100, 200),
    end = c(300, 400), strand = c("+", "+"),
    id = c("same_id", "same_id"),
    stringsAsFactors = FALSE
  )
  write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)

  spec <- orf_caller_spec(
    name = "test_dup", file_format = "tsv",
    column_map = list(chrom = "chrom", start = "start", end = "end",
                      strand = "strand", orf_id = "id")
  )
  register_orf_caller(spec, overwrite = TRUE)

  gr <- parse_orfs(file, source = "test_dup")
  expect_false(any(duplicated(names(gr))))

  unlink(file)
})


# ==============================================================================
# 10. Integration with package data
# ==============================================================================

test_that("parse_orfs works on bundled GENCODE BED file", {
  bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
  skip_if(bed == "", message = "Package extdata not available")

  gr <- parse_orfs(bed, source = "gencode")
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) > 0)
  expect_true("orf_id" %in% colnames(mcols(gr)))
  expect_true("orf_type" %in% colnames(mcols(gr)))
})
