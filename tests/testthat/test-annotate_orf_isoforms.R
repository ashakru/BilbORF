test_that("annotate_orf_isoforms correctly annotates ORFs", {
  # Setup: Load test data
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
  skip_if(gtf == "" || bed == "", "Test data not found")
  
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(rtracklayer))
  BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  
  # Prepare annotations
  annotations <- suppressWarnings(suppressMessages(
    prepare_annotations_fromGTF(gtf, BSgenome)
  ))
  
  # Load ORFs
  orfs <- import(bed, format = "BED")
  names(orfs) <- orfs$name
  
  # Prepare transcript metadata
  transcripts_meta <- import(gtf, format = "GTF") %>%
    as.data.frame() %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::select(gene_name, gene_id, transcript_id, transcript_type) %>%
    dplyr::distinct()
  
  # Test function execution
  expect_no_error({
    annotated_orfs <- annotate_orf_isoforms(
      annotations, 
      orfs, 
      BSgenome,
      transcripts_meta
    )
  })
  
  # Run annotation
  annotated_orfs <- annotate_orf_isoforms(
    annotations, 
    orfs, 
    BSgenome,
    transcripts_meta
  )
  
  # Test output structure
  expect_type(annotated_orfs, "list")
  expect_named(annotated_orfs, c("table", "ranges"))
  
  # Test table structure
  expect_s3_class(annotated_orfs$table, "data.frame")
  expect_true(nrow(annotated_orfs$table) > 0)
  
  # Test required columns
  required_cols <- c("ORF_isoform_id", "ORF_id", "transcript_id", "gene_id",
                     "gene_name", "seq_nt", "seq_aa", "len_nt", "len_aa",
                     "start_codon", "stop_codon", "orf_status", "complete_codons")
  expect_true(all(required_cols %in% colnames(annotated_orfs$table)))
  
  # Test ranges structure
  expect_s4_class(annotated_orfs$ranges, "GRangesList")
  expect_true(length(annotated_orfs$ranges) > 0)
  expect_true(all(names(annotated_orfs$ranges) %in% annotated_orfs$table$ORF_isoform_id))
  
  # Test ORF status categories
  valid_statuses <- c("translatable", "internal_stop", "no_stop",
                      "no_start", "no_stop_no_start", "intronic", "no_transcript")
  expect_true(all(annotated_orfs$table$orf_status %in% valid_statuses))

  # Every input ORF must appear in the output at least once
  expect_true(all(names(orfs) %in% annotated_orfs$table$ORF_id))
  
  # Test that translatable ORFs have valid start/stop codons
  translatable <- annotated_orfs$table[annotated_orfs$table$orf_status == "translatable", ]
  if (nrow(translatable) > 0) {
    default_starts <- c("ATG", "TTG", "CTG", "GTG")
    default_stops <- c("TAG", "TAA", "TGA")
    expect_true(all(translatable$start_codon %in% default_starts))
    expect_true(all(translatable$stop_codon %in% default_stops))
  }
})

test_that("annotate_orf_isoforms handles custom start/stop codons", {
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
  
  # Test with custom codons
  annotated_orfs <- annotate_orf_isoforms(
    annotations, 
    orfs[1:5], # Use subset for speed
    BSgenome,
    transcripts_meta,
    start_codons = c("ATG"),  # Only canonical start
    stop_codons = c("TAA", "TAG", "TGA")
  )
  
  expect_s3_class(annotated_orfs$table, "data.frame")
  expect_true(nrow(annotated_orfs$table) > 0)
})

test_that("annotate_orf_isoforms validates input", {
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  skip_if(gtf == "", "Test data not found")
  
  annotations <- suppressWarnings(suppressMessages(
    prepare_annotations_fromGTF(gtf, BSgenome)
  ))
  
  transcripts_meta <- data.frame(
    transcript_id = names(annotations$transcripts)[1:5],
    gene_name = paste0("GENE", 1:5),
    gene_id = paste0("ENSG", 1:5)
  )
  
  # Test with NULL ORFs
  expect_error(
    annotate_orf_isoforms(annotations, NULL, BSgenome, transcripts_meta)
  )
  
  # Test with empty GRanges
  empty_orfs <- GenomicRanges::GRanges()
  expect_error(
    annotate_orf_isoforms(annotations, empty_orfs, BSgenome, transcripts_meta)
  )
})

test_that("annotate_orf_isoforms sequences are valid", {
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
    orfs[1:10], # subset for speed
    BSgenome,
    transcripts_meta
  )
  
  # Stub rows (intronic / no_transcript) have no sequence — restrict checks to
  # rows that have actual sequences
  has_seq <- !is.na(annotated_orfs$table$seq_nt)

  # Check nucleotide sequences contain only valid bases
  valid_nt <- c("A", "T", "C", "G", "N")
  nt_bases <- unique(unlist(strsplit(annotated_orfs$table$seq_nt[has_seq], "")))
  expect_true(all(nt_bases %in% valid_nt))

  # Check amino acid sequences contain only valid AAs
  valid_aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
                "P", "Q", "R", "S", "T", "V", "W", "Y", "*", "X")
  aa_chars <- unique(unlist(strsplit(annotated_orfs$table$seq_aa[has_seq], "")))
  expect_true(all(aa_chars %in% valid_aa))

  # Check length relationships for rows with sequence
  seq_rows <- annotated_orfs$table[has_seq, ]
  expect_true(all(seq_rows$len_nt == nchar(seq_rows$seq_nt)))
  expect_true(all(seq_rows$len_aa <= nchar(seq_rows$seq_aa)))

  # complete_codons must be logical; NA only for stub rows
  expect_type(annotated_orfs$table$complete_codons, "logical")
  expect_true(all(is.na(annotated_orfs$table$complete_codons[!has_seq])))
  expect_true(all(seq_rows$complete_codons == (seq_rows$len_nt %% 3 == 0)))

  # translatable ORFs must have complete codons
  translatable_rows <- annotated_orfs$table[
    !is.na(annotated_orfs$table$orf_status) &
    annotated_orfs$table$orf_status == "translatable", ]
  if (nrow(translatable_rows) > 0) {
    expect_true(all(translatable_rows$complete_codons))
  }

  # All input ORFs appear in the output
  expect_true(all(names(orfs[1:10]) %in% annotated_orfs$table$ORF_id))
})