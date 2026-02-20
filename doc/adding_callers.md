# Tutorial: Adding ORF Caller Support to BilbORF

This guide explains how to add support for new ORF prediction tools and how to parse their outputs using the BilbORF parser system.

## Table of Contents

1. [Parser System Overview](#parser-system-overview)
2. [Adding a New Caller](#adding-a-new-caller)
3. [Parsing ORF Predictions](#parsing-orf-predictions)
4. [Complete Examples](#complete-examples)
5. [Best Practices](#best-practices)

---

## Parser System Overview

BilbORF uses a **spec-based registry system** for parsing ORF predictions from different callers. Each caller has a specification that defines:

- **Column mappings**: How to map caller-specific columns to standardized names
- **Extra columns**: Additional columns to preserve from the original output
- **Read function**: How to load and parse the caller's output files
- **Post-processing**: Optional transformations after reading

This extensible design allows you to add new callers without modifying core parsing logic.

### Currently Supported Callers

- **RiboTISH** (.txt files)
- **RiboCode** (.txt files)
- **PRICE** (.bed files)
- **Ribotricer** (.tsv files)
- **ORFquant** (RData files with GRangesList objects)
- **GENCODE** (.gtf annotation files)
- **RiboTIE** (.csv files)

---

## Adding a New Caller

### Step 1: Understand the Spec Components

A caller spec consists of four main components:

```r
list(
  required_columns = list(...),  # Column name mappings
  extra_cols = character(),       # Additional columns to preserve
  read_fn = function(path) {...}, # Function to read the file
  post_process_fn = NULL          # Optional post-processing
)
```

### Step 2: Create Column Mappings

The `required_columns` list maps standard BilbORF column names to the caller's column names:

```r
required_columns = list(
  seqnames = "Chromosome",        # Chromosome/sequence name
  start = "Start",                # Start coordinate
  end = "Stop",                   # End coordinate
  strand = "Strand",              # Strand ("+", "-", or "*")
  transcript_id = "TranscriptID", # Transcript identifier
  orf_id = "ORF_ID",             # Unique ORF identifier
  orf_type = "ORF_Type"          # ORF category (CDS, uORF, etc.)
)
```

**Note**: If a column doesn't exist in the original output, you can create it in the `read_fn` or `post_process_fn`.

### Step 3: Define Extra Columns

List any additional columns from the original output that you want to preserve:

```r
extra_cols = c(
  "pvalue",
  "score",
  "ORF_length",
  "predicted_expression"
)
```

These columns will be available in `output = "full"` mode but stripped in `output = "unified"` mode.

### Step 4: Write the Read Function

The `read_fn` takes a file path and returns a data.frame:

```r
read_fn = function(path) {
  # Read the file
  df <- read.table(path, header = TRUE, sep = "\t", 
                   stringsAsFactors = FALSE, comment.char = "")
  
  # Create any missing required columns
  df$orf_id <- paste0("ORF_", seq_len(nrow(df)))
  
  # Clean/transform data as needed
  df$strand <- ifelse(df$strand == "1", "+", "-")
  
  return(df)
}
```

**Key requirements**:
- Must return a data.frame
- Must include all columns specified in `required_columns`
- Should handle missing files gracefully
- Should validate data types

### Step 5: Write Post-Processing Function (Optional)

The `post_process_fn` takes a GRanges object and returns a modified GRanges:

```r
post_process_fn = function(gr) {
  # Example: Extract orf_type from a complex column
  if ("orf_category" %in% colnames(mcols(gr))) {
    mcols(gr)$orf_type <- sapply(mcols(gr)$orf_category, 
                                  function(x) strsplit(x, "_")[[1]][1])
  }
  
  # Example: Add computed columns
  mcols(gr)$orf_length <- width(gr)
  
  return(gr)
}
```

Use this for:
- Extracting information from complex columns
- Computing derived values
- Standardizing orf_type categories

### Step 6: Register the Spec

The spec should be registered when the package loads. Add it to `R/parse_orfs.R` in the `.onLoad` function:

```r
register_orf_caller_spec(
  caller = "mycaller",
  spec = list(
    required_columns = list(
      seqnames = "chr",
      start = "start_pos",
      end = "end_pos",
      strand = "strand",
      transcript_id = "tx_id",
      orf_id = "orf_name",
      orf_type = "category"
    ),
    extra_cols = c("score", "pvalue"),
    read_fn = function(path) {
      df <- read.table(path, header = TRUE, sep = "\t")
      df$orf_id <- paste0(df$tx_id, "_", df$orf_name)
      return(df)
    },
    post_process_fn = NULL
  )
)
```

---

## Parsing ORF Predictions

### Basic Usage

```r
# Load the package
library(BilbORF)

# Parse ORF predictions
results <- parse_orfs(
  path = "path/to/predictions.txt",
  source = "mycaller",
  output = "unified"
)
```

### Output Modes

The `output` parameter controls how much information is returned:

#### 1. `output = "full"` (Default)

Returns a data.frame with **all original columns** plus standardized coordinate information.

```r
results <- parse_orfs("predictions.txt", source = "ribotish", output = "full")

# Returns data.frame with columns:
# - All columns from required_columns
# - All columns from extra_cols
# - Original column names preserved where available
```

**Use when**: You need complete access to all caller-specific information and metadata.

#### 2. `output = "unified"`

Returns a **GRanges object** with only core metadata columns:
- `orf_id`: Unique ORF identifier
- `transcript_id`: Parent transcript ID
- `orf_type`: ORF category (CDS, uORF, dORF, etc.)
- `orf_width`: Width of the ORF in the transcript

```r
results <- parse_orfs("predictions.txt", source = "ribotish", output = "unified")

# Returns GRanges with genomic coordinates and 4 metadata columns
# Optimized for performance and downstream analysis
# Each range represents one exon/segment of an ORF
```

**Use when**: You need standardized output for comparative analysis across different callers or for efficient downstream processing.

**Note**: This mode strips all caller-specific columns and focuses on essential standardized information.

#### 3. `output = "granges"`

Returns a **minimal GRanges object** with genomic coordinates and `orf_id` only.

```r
results <- parse_orfs("predictions.txt", source = "ribotish", output = "granges")

# Returns GRanges with coordinates and only orf_id metadata
# No transcript_id, orf_type, or other columns
```

**Use when**: You only need coordinate information (e.g., for overlap analysis or visualization) and want maximum memory efficiency.

### Working with Different File Formats

The parser handles various file formats automatically based on the caller spec:

```r
# Text files (tab-delimited)
parse_orfs("ribotish_output.txt", source = "ribotish")

# BED files
parse_orfs("price_output.bed", source = "price")

# CSV files
parse_orfs("ribotie_output.csv", source = "ribotie")

# RData files (R objects)
parse_orfs("orfquant_results.RData", source = "orfquant")

# GTF annotation files
parse_orfs("gencode.v35.annotation.gtf", source = "gencode")
```

### Batch Processing Multiple Files

```r
# Process multiple files from the same caller
files <- c("sample1_ribotish.txt", "sample2_ribotish.txt", "sample3_ribotish.txt")

results_list <- lapply(files, function(f) {
  parse_orfs(f, source = "ribotish", output = "unified")
})

# Combine into single GRanges
combined <- do.call(c, results_list)
```

---

## Complete Examples

### Example 1: Adding Support for RiboCode

RiboCode outputs tab-delimited files with these columns:
- `ORF_ID`, `ORF_type`, `Chromosome`, `StartPosition`, `StopPosition`, 
- `Strand`, `ORFscore`, `pvalue`, etc.

```r
register_orf_caller_spec(
  caller = "ribocode",
  spec = list(
    required_columns = list(
      seqnames = "Chromosome",
      start = "StartPosition",
      end = "StopPosition",
      strand = "Strand",
      transcript_id = "transcript_id",
      orf_id = "ORF_ID",
      orf_type = "ORF_type"
    ),
    extra_cols = c("ORFscore", "pvalue", "ORF_length"),
    read_fn = function(path) {
      df <- read.table(path, header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, comment.char = "")
      
      # Extract transcript_id from ORF_ID
      df$transcript_id <- sub("_\\d+$", "", df$ORF_ID)
      
      # RiboCode uses numeric strand encoding
      df$Strand <- ifelse(df$Strand == "1", "+", "-")
      
      return(df)
    },
    post_process_fn = NULL
  )
)
```

### Example 2: Adding Support for a Custom Caller

Suppose you have a custom ORF caller that outputs CSV files:

```csv
gene_id,tx_name,chr,tx_start,tx_end,direction,orf_class,confidence
ENSG001,ENST001,chr1,1000,1500,forward,canonical,0.95
ENSG001,ENST001,chr1,800,900,forward,upstream,0.87
```

```r
register_orf_caller_spec(
  caller = "customcaller",
  spec = list(
    required_columns = list(
      seqnames = "chr",
      start = "tx_start",
      end = "tx_end",
      strand = "strand",
      transcript_id = "tx_name",
      orf_id = "orf_id",
      orf_type = "orf_class"
    ),
    extra_cols = c("confidence", "gene_id"),
    read_fn = function(path) {
      df <- read.csv(path, stringsAsFactors = FALSE)
      
      # Convert direction to strand
      df$strand <- ifelse(df$direction == "forward", "+", "-")
      
      # Create unique ORF IDs
      df$orf_id <- paste(df$tx_name, df$orf_class, 
                         df$tx_start, df$tx_end, sep = "_")
      
      return(df)
    },
    post_process_fn = function(gr) {
      # Standardize orf_type categories
      type_map <- c(
        "canonical" = "CDS",
        "upstream" = "uORF",
        "downstream" = "dORF",
        "internal" = "intORF"
      )
      mcols(gr)$orf_type <- type_map[mcols(gr)$orf_type]
      return(gr)
    }
  )
)
```

### Example 3: Handling RData Files (like ORFquant)

```r
register_orf_caller_spec(
  caller = "orfquant",
  spec = list(
    required_columns = list(
      seqnames = "seqnames",
      start = "start",
      end = "end",
      strand = "strand",
      transcript_id = "transcript_id",
      orf_id = "orf_id",
      orf_type = "orf_type"
    ),
    extra_cols = c("ORF_category_Tx", "ORF_category_Tx_compatible", 
                   "ORF_category_Gen"),
    read_fn = function(path) {
      # Load RData file
      env <- new.env()
      load(path, envir = env)
      
      # ORFquant returns a list with ORFs_tx and ORFs_gen
      res <- env[[ls(env)[1]]]
      
      # Prefer transcript-based coordinates (has metadata)
      if ("ORFs_tx" %in% names(res)) {
        gr <- res$ORFs_tx
      } else if ("ORFs_gen" %in% names(res)) {
        gr <- res$ORFs_gen
      } else {
        stop("No recognized ORF data in ORFquant results")
      }
      
      # Convert GRangesList to data.frame
      df <- as.data.frame(gr)
      return(df)
    },
    post_process_fn = function(gr) {
      # Extract orf_type from category columns
      if ("ORF_category_Tx" %in% colnames(mcols(gr))) {
        mcols(gr)$orf_type <- mcols(gr)$ORF_category_Tx
      } else if ("ORF_category_Gen" %in% colnames(mcols(gr))) {
        mcols(gr)$orf_type <- mcols(gr)$ORF_category_Gen
      }
      return(gr)
    }
  )
)
```

---

## Best Practices

### 1. **Validate Input Data**

```r
read_fn = function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  
  df <- read.table(path, header = TRUE, sep = "\t")
  
  # Check for required columns before mapping
  required <- c("chr", "start", "end", "strand")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  return(df)
}
```

### 2. **Handle Edge Cases**

```r
read_fn = function(path) {
  df <- read.table(path, header = TRUE, sep = "\t")
  
  # Handle empty files
  if (nrow(df) == 0) {
    warning("Empty file: ", path)
    return(data.frame())
  }
  
  # Handle missing strand information
  if (!"strand" %in% colnames(df)) {
    warning("No strand information, setting to '*'")
    df$strand <- "*"
  }
  
  return(df)
}
```

### 3. **Standardize ORF Types**

Use consistent naming across callers:
- `CDS` - Canonical coding sequence
- `uORF` - Upstream ORF (5' UTR)
- `dORF` - Downstream ORF (3' UTR)
- `intORF` - Internal/overlapping ORF
- `ext` - Extension (N-terminal or C-terminal)
- `uoORF` - Upstream overlapping ORF

```r
post_process_fn = function(gr) {
  # Map caller-specific types to standard types
  type_map <- c(
    "annotated" = "CDS",
    "novel_uORF" = "uORF",
    "novel_dORF" = "dORF",
    "extended" = "ext"
  )
  
  mcols(gr)$orf_type <- type_map[as.character(mcols(gr)$orf_type)]
  
  # Handle unmapped types
  unmapped <- is.na(mcols(gr)$orf_type)
  if (any(unmapped)) {
    warning("Unmapped ORF types found, keeping original")
    mcols(gr)$orf_type[unmapped] <- mcols(gr)$original_type[unmapped]
  }
  
  return(gr)
}
```

### 4. **Create Unique ORF IDs**

Ensure each ORF has a unique identifier:

```r
read_fn = function(path) {
  df <- read.table(path, header = TRUE, sep = "\t")
  
  # Create unique IDs if not provided
  if (!"orf_id" %in% colnames(df)) {
    df$orf_id <- paste(
      df$transcript_id,
      df$orf_type,
      df$start,
      df$end,
      sep = "_"
    )
  }
  
  # Check for duplicates
  if (any(duplicated(df$orf_id))) {
    warning("Duplicate ORF IDs found, making unique")
    df$orf_id <- make.unique(df$orf_id)
  }
  
  return(df)
}
```

### 5. **Document Caller-Specific Quirks**

Add comments explaining any non-obvious transformations:

```r
read_fn = function(path) {
  df <- read.table(path, header = TRUE, sep = "\t")
  
  # NOTE: RiboCode uses 1-based coordinates (already compatible with GRanges)
  # No adjustment needed for start position
  
  # NOTE: RiboCode encodes strand as 1/0 instead of +/-
  df$Strand <- ifelse(df$Strand == "1", "+", "-")
  
  # NOTE: ORF_ID format is "ENST00000123456_1" where _1 is the ORF number
  # Extract transcript ID by removing the suffix
  df$transcript_id <- sub("_\\d+$", "", df$ORF_ID)
  
  return(df)
}
```

### 6. **Test with Real Data**

Always test your spec with actual output files:

```r
# Test parsing
test_file <- "path/to/real/output.txt"
result_full <- parse_orfs(test_file, source = "mycaller", output = "full")
result_unified <- parse_orfs(test_file, source = "mycaller", output = "unified")
result_granges <- parse_orfs(test_file, source = "mycaller", output = "granges")

# Verify structure
stopifnot(is.data.frame(result_full))
stopifnot(is(result_unified, "GRanges"))
stopifnot(is(result_granges, "GRanges"))

# Check required columns in unified output
required_cols <- c("orf_id", "transcript_id", "orf_type", "orf_width")
stopifnot(all(required_cols %in% colnames(mcols(result_unified))))
```

### 7. **Performance Considerations**

For large files or complex data structures:

```r
read_fn = function(path) {
  # Use data.table for large files
  if (file.size(path) > 1e8) {  # > 100MB
    df <- data.table::fread(path, data.table = FALSE)
  } else {
    df <- read.table(path, header = TRUE, sep = "\t")
  }
  
  return(df)
}

post_process_fn = function(gr) {
  # For unified output, strip unnecessary columns early
  # This happens automatically in the pipeline, but you can
  # do additional filtering here if needed
  
  return(gr)
}
```

---

## Summary

1. **Understand your caller's output format**: Column names, data types, coordinate system
2. **Create a spec**: Define column mappings, extra columns, read function, and post-processing
3. **Register the spec**: Add it to `.onLoad()` in `R/parse_orfs.R`
4. **Test thoroughly**: Use real data files and verify all output modes
5. **Document**: Add comments explaining any transformations or quirks

The spec-based system makes it easy to add new callers without modifying core code. Each caller's logic is self-contained and maintainable.

---

## Additional Resources

- See `R/parse_orfs.R` for complete implementations of all supported callers
- Check `tests/testthat/test-parse_orfs.R` for testing examples
- Review `inst/extdata/` for example test files

For questions or to contribute new caller support, please open an issue on the GitHub repository.
