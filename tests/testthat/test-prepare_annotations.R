test_that("prepare_annotations_fromGTF creates valid annotations", {
  # Skip if GTF file not available
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  skip_if(gtf == "", "GTF file not found")
  
  # Load BSgenome
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  
  # Test function execution
  expect_no_error({
    annotations <- prepare_annotations_fromGTF(gtf, BSgenome)
  })
  
  # Test output structure
  annotations <- suppressWarnings(suppressMessages(
    prepare_annotations_fromGTF(gtf, BSgenome)
  ))
  
  expect_type(annotations, "list")
  expect_named(annotations, "transcripts")
  expect_s4_class(annotations$transcripts, "GRangesList")
  
  # Test that transcripts have names
  expect_true(length(names(annotations$transcripts)) > 0)
  expect_true(all(nchar(names(annotations$transcripts)) > 0))
  
  # Test that transcripts contain genomic ranges
  expect_s4_class(annotations$transcripts[[1]], "GRanges")
  expect_true(length(annotations$transcripts[[1]]) > 0)
})

test_that("prepare_annotations_fromGTF handles invalid inputs", {
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  
  # Test with non-existent file
  expect_error(
    prepare_annotations_fromGTF("nonexistent.gtf", BSgenome)
  )
  
  # Test with NULL BSgenome
  gtf <- system.file("extdata", "gencode.v35.annotation_chr10.gtf", package = "BilbORF")
  skip_if(gtf == "", "GTF file not found")
  
  expect_error(
    prepare_annotations_fromGTF(gtf, NULL)
  )
})
