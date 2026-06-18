# ==============================================================================
# ORF Caller Spec System
#
# Extensible architecture for parsing ORF caller outputs. Each caller is
# described by an `orf_caller_spec` S3 object that encodes file format,
# column mappings, coordinate conventions, and optional custom logic.
#
# Users can add support for new callers via `register_orf_caller()` without
# modifying any existing code.
# ==============================================================================


# -- Registry environment (package-private) ------------------------------------

#' @keywords internal
.orf_caller_registry <- new.env(parent = emptyenv())


# -- S3 constructor & validator ------------------------------------------------

#' Create an ORF Caller Specification
#'
#' Defines how to read and interpret the output of a specific ORF caller.
#' The spec is used by \code{\link{parse_orfs}} to convert caller-specific
#' files into a standardised GRanges object.
#'
#' @param name Short identifier (lowercase, no spaces). Used as the
#'   \code{source} argument in \code{parse_orfs()}.
#' @param description One-line human-readable description of the caller.
#' @param file_format One of \code{"tsv"}, \code{"csv"}, \code{"bed"},
#'   \code{"rds"}. Controls how the file is initially loaded.
#' @param column_map A named list that maps **canonical field names** to the
#'   caller-specific column name(s). The canonical fields are:
#'   \describe{
#'     \item{chrom}{(required) Chromosome / seqname column.}
#'     \item{start}{(required) Start coordinate column.}
#'     \item{end}{(required) End coordinate column.}
#'     \item{strand}{(required) Strand column.}
#'     \item{orf_id}{(required) Unique ORF identifier column.}
#'     \item{orf_type}{(optional) ORF category column (e.g. uORF, dORF, CDS).}
#'     \item{gene_id}{(optional) Gene identifier column.}
#'     \item{gene_name}{(optional) Gene symbol column.}
#'     \item{transcript_id}{(optional) Transcript identifier column.}
#'   }
#'   Each value can be either:
#'   \itemize{
#'     \item A single string: The column name in the caller output.
#'     \item A character vector: Column name aliases tried in order (first
#'       match wins).
#'   }
#' @param coord_system Coordinate convention of the caller output.
#'   \code{"1-based"} (default, used by R/GRanges) or \code{"0-based"}
#'   (BED convention; start is adjusted +1 internally).
#' @param read_fn Optional custom function \code{function(file) -> data.frame}
#'   that overrides the default file reader. Useful when the file needs
#'   special parsing (e.g. composite columns like Ribo-TISH's GenomePos).
#'   Must return a data.frame, GRanges, or GRangesList.
#' @param post_process_fn Optional function
#'   \code{function(gr, raw_df) -> GRanges} called after the standard
#'   pipeline to add caller-specific metadata or transformations.
#'   Receives the GRanges built so far and the raw data.frame.
#' @param url URL to the caller's homepage or documentation.
#' @param extra_cols Character vector of additional caller-specific columns
#'   worth retaining by default (beyond the canonical fields).
#'
#' @return An S3 object of class \code{orf_caller_spec}.
#'
#' @export
#'
#' @examples
#' # Minimal spec for a hypothetical new caller
#' my_spec <- orf_caller_spec(
#'   name        = "my_caller",
#'   description = "My custom ORF caller",
#'   file_format = "tsv",
#'   column_map  = list(
#'     chrom  = "chr",
#'     start  = "orf_start",
#'     end    = "orf_end",
#'     strand = "strand",
#'     orf_id = "id"
#'   )
#' )
#'
#' # Register it so parse_orfs() can use it
#' register_orf_caller(my_spec)
#' print(supported_orf_callers())
#'
#' @seealso \code{\link{register_orf_caller}}, \code{\link{parse_orfs}}
orf_caller_spec <- function(name,
                            description     = "",
                            file_format     = c("tsv", "csv", "bed", "rds"),
                            column_map      = list(),
                            coord_system    = c("1-based", "0-based"),
                            read_fn         = NULL,
                            post_process_fn = NULL,
                            url             = "",
                            extra_cols      = character(0)) {

  file_format  <- match.arg(file_format)
  coord_system <- match.arg(coord_system)

  # Validate name
  if (!is.character(name) || length(name) != 1 || nchar(name) == 0) {
    stop("`name` must be a non-empty string.")
  }

  # Validate column_map -- require coordinate fields for tabular formats
  # unless a custom read_fn handles everything
  if (file_format %in% c("tsv", "csv") && is.null(read_fn)) {
    required_fields <- c("chrom", "start", "end", "strand", "orf_id")
    missing <- setdiff(required_fields, names(column_map))
    if (length(missing) > 0) {
      stop("column_map is missing required canonical fields: ",
           paste(missing, collapse = ", "),
           ". Required: ", paste(required_fields, collapse = ", "))
    }
  }

  # Validate callables
  if (!is.null(read_fn) && !is.function(read_fn)) {
    stop("`read_fn` must be a function or NULL.")
  }
  if (!is.null(post_process_fn) && !is.function(post_process_fn)) {
    stop("`post_process_fn` must be a function or NULL.")
  }

  spec <- structure(
    list(
      name            = tolower(name),
      description     = description,
      file_format     = file_format,
      column_map      = column_map,
      coord_system    = coord_system,
      read_fn         = read_fn,
      post_process_fn = post_process_fn,
      url             = url,
      extra_cols      = extra_cols
    ),
    class = "orf_caller_spec"
  )

  spec
}


#' @export
print.orf_caller_spec <- function(x, ...) {
  cat(sprintf("ORF Caller Spec: %s\n", x$name))
  cat(sprintf("  Description : %s\n", x$description))
  cat(sprintf("  File format : %s\n", x$file_format))
  cat(sprintf("  Coordinates : %s\n", x$coord_system))
  mapped <- names(x$column_map)
  cat(sprintf("  Mapped fields: %s\n", paste(mapped, collapse = ", ")))
  if (nchar(x$url) > 0) cat(sprintf("  URL         : %s\n", x$url))
  if (!is.null(x$read_fn)) cat("  Custom reader : yes\n")
  if (!is.null(x$post_process_fn)) cat("  Post-process  : yes\n")
  invisible(x)
}


# -- Registry functions --------------------------------------------------------

#' Register an ORF Caller
#'
#' Adds an \code{\link{orf_caller_spec}} to the internal registry so it can
#' be used with \code{\link{parse_orfs}}.
#'
#' @param spec An \code{orf_caller_spec} object.
#' @param overwrite Logical. If \code{TRUE}, silently replace an existing
#'   spec with the same name. Default \code{FALSE} (warns on overwrite).
#'
#' @return Invisibly returns the registered spec.
#' @export
#'
#' @examples
#' spec <- orf_caller_spec(
#'   name        = "my_caller",
#'   description = "My custom ORF caller",
#'   file_format = "tsv",
#'   column_map  = list(
#'     chrom  = "chr",
#'     start  = "orf_start",
#'     end    = "orf_end",
#'     strand = "strand",
#'     orf_id = "id"
#'   )
#' )
#' register_orf_caller(spec)
#'
#' @seealso \code{\link{orf_caller_spec}}, \code{\link{get_orf_caller}},
#'   \code{\link{supported_orf_callers}}
register_orf_caller <- function(spec, overwrite = FALSE) {
  if (!inherits(spec, "orf_caller_spec")) {
    stop("`spec` must be an orf_caller_spec object.")
  }

  if (exists(spec$name, envir = .orf_caller_registry) && !overwrite) {
    warning("Overwriting existing spec for '", spec$name,
            "'. Use overwrite = TRUE to suppress this warning.")
  }

  assign(spec$name, spec, envir = .orf_caller_registry)
  invisible(spec)
}


#' Get a Registered ORF Caller Spec
#'
#' Retrieves an \code{orf_caller_spec} from the registry by name.
#'
#' @param name Character string. The caller name (case-insensitive).
#' @return An \code{orf_caller_spec} object, or an error if not found.
#' @export
#'
#' @examples
#' # List available callers first
#' supported_orf_callers()
#'
#' @seealso \code{\link{register_orf_caller}}, \code{\link{supported_orf_callers}}
get_orf_caller <- function(name) {
  name <- tolower(name)
  if (!exists(name, envir = .orf_caller_registry)) {
    available <- ls(.orf_caller_registry)
    stop("Unknown ORF caller '", name, "'. ",
         "Available: ", paste(available, collapse = ", "), ". ",
         "Use register_orf_caller() to add new callers.")
  }
  get(name, envir = .orf_caller_registry)
}


#' List Supported ORF Callers
#'
#' Returns a data.frame describing all registered ORF callers.
#'
#' @return A data.frame with columns: source, description, file_format,
#'   coord_system, url.
#' @export
#'
#' @examples
#' supported_orf_callers()
supported_orf_callers <- function() {
  specs <- ls(.orf_caller_registry)
  if (length(specs) == 0) {
    return(data.frame(
      source = character(0), description = character(0),
      file_format = character(0), coord_system = character(0),
      url = character(0), stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(specs, function(nm) {
    s <- get(nm, envir = .orf_caller_registry)
    data.frame(
      source       = s$name,
      description  = s$description,
      file_format  = s$file_format,
      coord_system = s$coord_system,
      url          = s$url,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


# -- Generic parse engine ------------------------------------------------------

#' Parse ORF Caller Results into a Standardised GRanges Object
#'
#' A unified interface for importing ORF predictions from Ribo-seq ORF
#' callers. Converts caller-specific output formats into a named GRanges
#' object ready for use with \code{\link{annotate_orf_isoforms}}.
#'
#' New callers can be supported by creating an \code{\link{orf_caller_spec}}
#' and registering it with \code{\link{register_orf_caller}}.
#'
#' @param file Path to the ORF caller output file.
#' @param source Name of the ORF caller (case-insensitive). Must match a
#'   registered \code{orf_caller_spec}. See \code{supported_orf_callers()}.
#' @param genome_style \code{"UCSC"} (default, "chr1") or \code{"Ensembl"}
#'   ("1"). Chromosome names are converted to match.
#' @param min_length Minimum ORF length in nucleotides. ORFs shorter than
#'   this are dropped. Default \code{NULL} (no filter).
#' @param output Output format: \code{"unified"} (default) returns a
#'   GRanges with core metadata (genomic coordinates, transcript_id, orf_type,
#'   orf_width, orf_id); \code{"full"} returns a list with both
#'   the unified GRanges and a data.frame of all caller-specific metrics.
#' @param additional_cols Character vector of extra caller-specific columns
#'   to retain in the output. Columns are included in both \code{"unified"}
#'   and \code{"full"} output modes.
#' @param txdb Optional TxDb object or GRangesList of transcript exon structures.
#'   Required for callers like RiboCode that need splicing-aware coordinate
#'   adjustment. If NULL (default), genomic coordinate extension is used with
#'   a warning for multi-exon ORFs. Can be created with
#'   \code{GenomicFeatures::makeTxDbFromGFF()} or loaded from packages like
#'   \code{TxDb.Hsapiens.UCSC.hg38.knownGene}.
#' @param ribocode_extend_stop Logical. If \code{TRUE} (default), RiboCode
#'   coordinates are extended by 3 nt to include the stop codon (RiboCode reports
#'   coordinates excluding the stop codon). Set to \code{FALSE} to use
#'   RiboCode's original coordinates without modification. Ignored for other callers.
#'
#' @return If \code{output = "unified"} (default), a GRanges with one range
#'   per exon/segment and metadata columns: \code{orf_id} (unique identifier),
#'   \code{transcript_id}, \code{orf_type}, \code{orf_width} (width of this range).
#'   If \code{output = "full"}, a list with two elements: \code{orfs} (the unified
#'   GRanges) and \code{metadata} (data.frame with all caller metrics).
#'
#' @details
#' The parse pipeline:
#' \enumerate{
#'   \item Load the file using the spec's reader (\code{read_fn} or default
#'     for tsv/csv/bed/rds).
#'   \item Resolve column aliases via \code{column_map}.
#'   \item Build a GRanges using the canonical coordinate fields.
#'   \item Adjust coordinates if \code{coord_system = "0-based"} (+1 to start).
#'   \item Run \code{post_process_fn} if defined.
#'   \item Map additional metadata columns.
#'   \item Ensure unique names, apply genome style, filter by length.
#'   \item Format output: unified mode returns a GRanges with core metadata
#'     columns; full mode also returns all caller-specific metrics in a
#'     separate data.frame.
#' }
#'
#' The unified output is compatible with ORFik's suite of ORF analysis
#' functions for downstream processing.
#'
#' @export
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @importFrom GenomicRanges GRanges mcols mcols<-
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<- DataFrame
#'
#' @examples
#' # Parse GENCODE consensus ORFs (ships with package)
#' bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
#' if (bed != "") {
#'   # Unified output: GRangesList with core metadata
#'   orfs <- parse_orfs(bed, source = "gencode")
#'   print(length(orfs))  # Number of ORFs
#'   print(orfs[[1]])     # First ORF's ranges
#'   
#'   # Full output: includes all caller metrics
#'   full <- parse_orfs(bed, source = "gencode", output = "full")
#'   print(names(full))           # "orfs" and "metadata"
#'   print(head(full$metadata))   # All columns from caller
#' }
#'
#' \dontrun{
#' # Parse Ribo-TISH results (unified output by default)
#' orfs <- parse_orfs("ribotish_pred.txt", source = "ribotish")
#'
#' # Parse RiboCode with Ensembl chromosomes, get full metrics
#' result <- parse_orfs("ribocode.txt",
#'                      source = "ribocode",
#'                      genome_style = "Ensembl",
#'                      output = "full")
#' orfs_grl <- result$orfs
#' metrics <- result$metadata  # All RiboCode columns (pval, Psites, etc.)
#'
#' # Use ORFik functions on the unified output
#' library(ORFik)
#' orfs <- parse_orfs("ribocode.txt", source = "ribocode")
#' # Calculate total ORF widths per group
#' widths <- sum(width(orfs))
#' # Extract start codons (requires genome/transcript sequences)
#' # start_codons <- startCodons(orfs, faFile = "genome.fa")
#' 
#' # Register and use a custom caller
#' my_spec <- orf_caller_spec(
#'   name = "my_caller", file_format = "tsv",
#'   column_map = list(chrom = "chr", start = "begin",
#'                     end = "stop", strand = "str",
#'                     orf_id = "name")
#' )
#' register_orf_caller(my_spec)
#' orfs <- parse_orfs("my_results.tsv", source = "my_caller")
#' }
#'
#' @seealso \code{\link{orf_caller_spec}}, \code{\link{register_orf_caller}},
#'   \code{\link{supported_orf_callers}},
#'   \code{\link{annotate_orf_isoforms}}
parse_orfs <- function(file,
                       source,
                       genome_style    = "UCSC",
                       min_length      = NULL,
                       output          = c("unified", "full"),
                       additional_cols = NULL,
                       txdb            = NULL,
                       ribocode_extend_stop = TRUE) {

  # --- Validate inputs --------------------------------------------------------
  if (!file.exists(file)) stop("File not found: ", file)
  if (!genome_style %in% c("UCSC", "Ensembl")) {
    stop("genome_style must be 'UCSC' or 'Ensembl'.")
  }
  output <- match.arg(output)

  spec <- get_orf_caller(source)

  # --- 1. Read the file -------------------------------------------------------
  raw <- .read_file(file, spec)

  # --- 2-3. Build GRanges from raw data ---------------------------------------
  gr <- .build_granges(raw, spec, additional_cols)

  # --- 4. Post-process --------------------------------------------------------
  if (!is.null(spec$post_process_fn)) {
    # Post-process functions may need txdb for splicing-aware operations
    if (is.data.frame(raw)) {
      gr <- spec$post_process_fn(gr, raw, txdb, ribocode_extend_stop)
    } else {
      gr <- spec$post_process_fn(gr, NULL, txdb, ribocode_extend_stop)
    }
  }

  # --- 4.5. Strip metadata early for unified output (performance optimization) --
  # For unified output, we only need transcript_id, orf_type, orf_id (+ any
  # additional_cols requested by the user).
  # Stripping unnecessary columns (especially complex types like CompressedLists)
  # before genome style conversion significantly improves performance
  if (output == "unified") {
    cols_to_keep <- c("transcript_id", "orf_type", "orf_id", additional_cols)
    existing_cols <- cols_to_keep[cols_to_keep %in% colnames(GenomicRanges::mcols(gr))]
    if (length(existing_cols) > 0 && length(existing_cols) < ncol(GenomicRanges::mcols(gr))) {
      GenomicRanges::mcols(gr) <- GenomicRanges::mcols(gr)[, existing_cols, drop = FALSE]
    }
  }

  # --- 5. Ensure unique names --------------------------------------------------
  if (is.null(names(gr)) || any(names(gr) == "") || any(is.na(names(gr)))) {
    names(gr) <- paste0(spec$name, "_orf_", seq_along(gr))
  }
  if (any(duplicated(names(gr)))) {
    names(gr) <- make.unique(names(gr), sep = "_")
  }
  if (!"orf_id" %in% colnames(GenomicRanges::mcols(gr))) {
    GenomicRanges::mcols(gr)$orf_id <- names(gr)
  }

  # --- 6. Genome style --------------------------------------------------------
  tryCatch({
    GenomeInfoDb::seqlevelsStyle(gr) <- genome_style
  }, error = function(e) {
    warning("Could not set seqlevels style to '", genome_style, "': ",
            conditionMessage(e), ". Leaving names unchanged.")
  })

  # --- 7. Length filter --------------------------------------------------------
  if (!is.null(min_length)) {
    orf_len <- if (is(gr, "GRangesList")) sum(width(gr)) else width(gr)
    gr <- gr[orf_len >= min_length]
    if (length(gr) == 0) {
      warning("All ORFs filtered by min_length = ", min_length)
    }
  }

  # --- 8. Format output -------------------------------------------------------
  if (output == "unified") {
    return(.make_unified_output(gr, additional_cols))
  } else {
    return(.make_full_output(gr, raw))
  }
}


# ==============================================================================
# Combine multiple ORF datasets
# ==============================================================================

#' Collapse ORF Calls from Multiple Datasets
#'
#' Takes a list of parsed ORF caller outputs (unified mode) and creates a
#' combined GRanges with metadata columns indicating which datasets contain
#' each ORF. This is useful for comparing ORF predictions across multiple
#' samples, conditions, or replicates.
#'
#' @param orf_list A list of GRanges objects, typically from
#'   \code{\link{parse_orfs}} with \code{output = "unified"}.
#' @param dataset_names Optional character vector of names for each dataset.
#'   If NULL, uses list names or generates "dataset_1", "dataset_2", etc.
#' @param match_by How to determine if ORFs are "the same" across datasets:
#'   \describe{
#'     \item{coordinates}{ORFs must have identical genomic coordinates (default)}
#'     \item{start_codon}{ORFs match if they have the same start position
#'                        (seqnames, start, strand), ignoring transcript ID}
#'     \item{start_codon_transcript}{ORFs match if they have the same start
#'                                   position and transcript_id}
#'     \item{stop_codon}{ORFs match if they have the same stop position
#'                       (seqnames, end, strand), ignoring transcript ID}
#'     \item{stop_codon_transcript}{ORFs match if they have the same stop
#'                                  position and transcript_id}
#'     \item{start_stop}{ORFs match if they have the same start and stop
#'                       positions (seqnames, start, end, strand)}
#'     \item{start_stop_transcript}{ORFs match if they have the same start, stop,
#'                                  and transcript_id}
#'     \item{overlap}{ORFs must overlap by at least \code{overlap_threshold}}
#'     \item{metadata}{ORFs match if they have the same values for
#'                     \code{metadata_keys} (e.g., same transcript + orf_type)}
#'   }
#' @param overlap_threshold Numeric in [0, 1]. For \code{match_by = "overlap"},
#'   the minimum fraction of reciprocal overlap required (default: 0.9).
#' @param metadata_keys Character vector of metadata column names to use for
#'   matching when \code{match_by = "metadata"}. Default: c("transcript_id", "orf_type").
#' @param combine How to combine ORFs across datasets:
#'   \describe{
#'     \item{union}{Return all unique ORFs from all datasets (default)}
#'     \item{first}{Use first dataset as reference, add presence indicators}
#'   }
#' @param keep_orf_ids Logical. If TRUE, retain the ORF names and orf_id metadata
#'   column from the input GRanges. If FALSE (default), drop both ORF names and
#'   the orf_id column from the result, since combined ORFs from different datasets
#'   may have conflicting IDs.
#'
#' @return A GRanges object with:
#'   \itemize{
#'     \item All metadata columns from the input GRanges
#'     \item Additional logical columns named \code{in_<dataset_name>} indicating
#'           presence in each dataset
#'     \item \code{n_datasets}: Integer column with count of datasets containing
#'           this ORF
#'   }
#'
#' @examples
#' \dontrun{
#' # Parse ORFquant results from multiple samples
#' orfquant_files <- c("sample1.rds", "sample2.rds", "sample3.rds")
#' orfquant_orfs <- lapply(orfquant_files, function(f) {
#'   parse_orfs(f, source = "orfquant", genome_style = "UCSC")
#' })
#' names(orfquant_orfs) <- c("control", "treatment_1", "treatment_2")
#' 
#' # Collapse by exact coordinates
#' combined <- collapse_orf_calls(orfquant_orfs)
#' 
#' # See which ORFs are in all datasets
#' table(combined$n_datasets)
#' common_orfs <- combined[combined$n_datasets == 3]
#' 
#' # Keep original ORF IDs from input datasets
#' combined_with_ids <- collapse_orf_calls(orfquant_orfs, keep_orf_ids = TRUE)
#' head(names(combined_with_ids))  # Shows ORF IDs
#' 
#' # Match by start codon only (ignores transcript ID and stop differences)
#' combined_start <- collapse_orf_calls(
#'   orfquant_orfs,
#'   match_by = "start_codon"
#' )
#' 
#' # Match by start codon and transcript (allows stop codon differences)
#' combined_start_tx <- collapse_orf_calls(
#'   orfquant_orfs,
#'   match_by = "start_codon_transcript"
#' )
#' 
#' # Match by stop codon only
#' combined_stop <- collapse_orf_calls(
#'   orfquant_orfs,
#'   match_by = "stop_codon"
#' )
#' 
#' # Match by stop codon and transcript
#' combined_stop_tx <- collapse_orf_calls(
#'   orfquant_orfs,
#'   match_by = "stop_codon_transcript"
#' )
#' 
#' # Match by start and stop positions
#' combined_start_stop <- collapse_orf_calls(
#'   orfquant_orfs,
#'   match_by = "start_stop"
#' )
#' 
#' # Match by start, stop, and transcript ID
#' combined_exact <- collapse_orf_calls(
#'   orfquant_orfs,
#'   match_by = "start_stop_transcript"
#' )
#' 
#' # Collapse by overlap (allows slight coordinate differences)
#' combined_overlap <- collapse_orf_calls(
#'   orfquant_orfs,
#'   match_by = "overlap",
#'   overlap_threshold = 0.95
#' )
#' 
#' # Collapse by metadata (same transcript + orf_type)
#' combined_meta <- collapse_orf_calls(
#'   orfquant_orfs,
#'   match_by = "metadata",
#'   metadata_keys = c("transcript_id", "orf_type")
#' )
#' }
#'
#' @export
collapse_orf_calls <- function(orf_list,
                                dataset_names = NULL,
                                match_by = c("coordinates", "start_codon", "start_codon_transcript",
                                            "stop_codon", "stop_codon_transcript",
                                            "start_stop", "start_stop_transcript",
                                            "overlap", "metadata"),
                                overlap_threshold = 0.9,
                                metadata_keys = c("transcript_id", "orf_type"),
                                combine = c("union", "first"),
                                keep_orf_ids = FALSE) {
  
  # --- 1. Input validation ----------------------------------------------------
  if (!is.list(orf_list) || length(orf_list) < 1) {
    stop("orf_list must be a list with at least one element")
  }
  
  # Check all elements are GRanges
  if (!all(vapply(orf_list, function(x) is(x, "GRanges"), logical(1)))) {
    stop("All elements of orf_list must be GRanges objects")
  }
  
  match_by <- match.arg(match_by)
  combine <- match.arg(combine)
  
  # Handle dataset names (before early return for single dataset)
  n_datasets <- length(orf_list)
  if (is.null(dataset_names)) {
    if (!is.null(names(orf_list))) {
      dataset_names <- names(orf_list)
    } else {
      dataset_names <- paste0("dataset_", seq_len(n_datasets))
    }
  } else {
    if (length(dataset_names) != n_datasets) {
      stop("dataset_names must have the same length as orf_list")
    }
  }
  
  # Ensure valid R column names
  dataset_names <- make.names(dataset_names, unique = TRUE)
  
  # Handle single dataset case
  if (length(orf_list) == 1) {
    warning("Only one dataset provided. Returning it with n_datasets = 1")
    gr <- orf_list[[1]]
    GenomicRanges::mcols(gr)$n_datasets <- 1L
    col_name <- paste0("in_", dataset_names[1])
    GenomicRanges::mcols(gr)[[col_name]] <- TRUE
    # Drop ORF IDs if requested
    if (!keep_orf_ids) {
      names(gr) <- NULL
      # Also remove orf_id column from metadata
      if ("orf_id" %in% colnames(GenomicRanges::mcols(gr))) {
        GenomicRanges::mcols(gr)$orf_id <- NULL
      }
    }
    return(gr)
  }
  
  # --- 2. Get base set of ORFs ------------------------------------------------
  if (combine == "first") {
    base_gr <- orf_list[[1]]
  } else {  # union
    # Combine all, keeping duplicates for now
    base_gr <- do.call(c, unname(orf_list))
  }
  
  # --- 3. Match ORFs across datasets ------------------------------------------
  
  # Initialize presence matrix
  presence_matrix <- matrix(FALSE, nrow = length(base_gr), ncol = n_datasets)
  colnames(presence_matrix) <- paste0("in_", dataset_names)
  
  if (match_by == "coordinates") {
    # Exact coordinate matching
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      # Find exact matches
      matches <- GenomicRanges::match(base_gr, query_gr)
      presence_matrix[, i] <- !is.na(matches)
    }
    
  } else if (match_by == "start_codon") {
    # Match by start position only (seqnames + start + strand)
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      
      # Create keys for start position
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::start(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        sep = "|"
      )
      
      query_key <- paste(
        as.character(GenomicRanges::seqnames(query_gr)),
        GenomicRanges::start(query_gr),
        as.character(GenomicRanges::strand(query_gr)),
        sep = "|"
      )
      
      # Match by key
      matches <- match(base_key, query_key)
      presence_matrix[, i] <- !is.na(matches)
    }
    
  } else if (match_by == "start_codon_transcript") {
    # Match by start position and transcript_id
    # Check that transcript_id exists
    if (!all(vapply(orf_list, function(gr) {
      "transcript_id" %in% colnames(GenomicRanges::mcols(gr))
    }, logical(1)))) {
      stop("transcript_id not found in all datasets. Required for match_by='start_codon_transcript'")
    }
    
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      
      # Create keys for start position and transcript
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::start(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        as.character(GenomicRanges::mcols(base_gr)$transcript_id),
        sep = "|"
      )
      
      query_key <- paste(
        as.character(GenomicRanges::seqnames(query_gr)),
        GenomicRanges::start(query_gr),
        as.character(GenomicRanges::strand(query_gr)),
        as.character(GenomicRanges::mcols(query_gr)$transcript_id),
        sep = "|"
      )
      
      # Match by key
      matches <- match(base_key, query_key)
      presence_matrix[, i] <- !is.na(matches)
    }
    
  } else if (match_by == "stop_codon") {
    # Match by stop position only (seqnames + end + strand)
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      
      # Create keys for stop position
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::end(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        sep = "|"
      )
      
      query_key <- paste(
        as.character(GenomicRanges::seqnames(query_gr)),
        GenomicRanges::end(query_gr),
        as.character(GenomicRanges::strand(query_gr)),
        sep = "|"
      )
      
      # Match by key
      matches <- match(base_key, query_key)
      presence_matrix[, i] <- !is.na(matches)
    }
    
  } else if (match_by == "stop_codon_transcript") {
    # Match by stop position and transcript_id
    # Check that transcript_id exists
    if (!all(vapply(orf_list, function(gr) {
      "transcript_id" %in% colnames(GenomicRanges::mcols(gr))
    }, logical(1)))) {
      stop("transcript_id not found in all datasets. Required for match_by='stop_codon_transcript'")
    }
    
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      
      # Create keys for stop position and transcript
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::end(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        as.character(GenomicRanges::mcols(base_gr)$transcript_id),
        sep = "|"
      )
      
      query_key <- paste(
        as.character(GenomicRanges::seqnames(query_gr)),
        GenomicRanges::end(query_gr),
        as.character(GenomicRanges::strand(query_gr)),
        as.character(GenomicRanges::mcols(query_gr)$transcript_id),
        sep = "|"
      )
      
      # Match by key
      matches <- match(base_key, query_key)
      presence_matrix[, i] <- !is.na(matches)
    }
    
  } else if (match_by == "start_stop") {
    # Match by start and stop positions (seqnames + start + end + strand)
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      
      # Create keys for start and stop positions
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::start(base_gr),
        GenomicRanges::end(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        sep = "|"
      )
      
      query_key <- paste(
        as.character(GenomicRanges::seqnames(query_gr)),
        GenomicRanges::start(query_gr),
        GenomicRanges::end(query_gr),
        as.character(GenomicRanges::strand(query_gr)),
        sep = "|"
      )
      
      # Match by key
      matches <- match(base_key, query_key)
      presence_matrix[, i] <- !is.na(matches)
    }
    
  } else if (match_by == "start_stop_transcript") {
    # Match by start, stop, and transcript_id
    # Check that transcript_id exists
    if (!all(vapply(orf_list, function(gr) {
      "transcript_id" %in% colnames(GenomicRanges::mcols(gr))
    }, logical(1)))) {
      stop("transcript_id not found in all datasets. Required for match_by='start_stop_transcript'")
    }
    
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      
      # Create keys for start, stop, and transcript
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::start(base_gr),
        GenomicRanges::end(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        as.character(GenomicRanges::mcols(base_gr)$transcript_id),
        sep = "|"
      )
      
      query_key <- paste(
        as.character(GenomicRanges::seqnames(query_gr)),
        GenomicRanges::start(query_gr),
        GenomicRanges::end(query_gr),
        as.character(GenomicRanges::strand(query_gr)),
        as.character(GenomicRanges::mcols(query_gr)$transcript_id),
        sep = "|"
      )
      
      # Match by key
      matches <- match(base_key, query_key)
      presence_matrix[, i] <- !is.na(matches)
    }
    
  } else if (match_by == "overlap") {
    # Overlap-based matching
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      # Find overlaps
      hits <- GenomicRanges::findOverlaps(base_gr, query_gr)
      
      # Calculate reciprocal overlap
      query_idx <- S4Vectors::subjectHits(hits)
      base_idx <- S4Vectors::queryHits(hits)
      
      overlap_width <- IRanges::width(GenomicRanges::pintersect(
        base_gr[base_idx], query_gr[query_idx]
      ))
      base_width <- IRanges::width(base_gr[base_idx])
      query_width <- IRanges::width(query_gr[query_idx])
      
      # Reciprocal overlap fraction
      overlap_frac <- pmin(overlap_width / base_width,
                          overlap_width / query_width)
      
      # Mark as present if overlap threshold met
      good_hits <- base_idx[overlap_frac >= overlap_threshold]
      presence_matrix[unique(good_hits), i] <- TRUE
    }
    
  } else {  # metadata
    # Check metadata columns exist
    for (key in metadata_keys) {
      if (!all(vapply(orf_list, function(gr) {
        key %in% colnames(GenomicRanges::mcols(gr))
      }, logical(1)))) {
        stop("metadata_key '", key, "' not found in all datasets")
      }
    }
    
    # Build metadata keys for each dataset
    for (i in seq_len(n_datasets)) {
      query_gr <- orf_list[[i]]
      
      # Create composite key for base and query
      base_key <- do.call(paste, c(
        lapply(metadata_keys, function(k) {
          as.character(GenomicRanges::mcols(base_gr)[[k]])
        }),
        sep = "|"
      ))
      
      query_key <- do.call(paste, c(
        lapply(metadata_keys, function(k) {
          as.character(GenomicRanges::mcols(query_gr)[[k]])
        }),
        sep = "|"
      ))
      
      # Match by key
      matches <- match(base_key, query_key)
      presence_matrix[, i] <- !is.na(matches)
    }
  }
  
  # --- 4. Remove duplicates if combine = "union" ------------------------------
  if (combine == "union") {
    # Find unique ORFs based on match_by strategy
    if (match_by == "coordinates") {
      is_unique <- !duplicated(base_gr)
      
    } else if (match_by == "start_codon") {
      # Unique by start position
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::start(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        sep = "|"
      )
      is_unique <- !duplicated(base_key)
      
    } else if (match_by == "start_codon_transcript") {
      # Unique by start position and transcript
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::start(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        as.character(GenomicRanges::mcols(base_gr)$transcript_id),
        sep = "|"
      )
      is_unique <- !duplicated(base_key)
      
    } else if (match_by == "stop_codon") {
      # Unique by stop position
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::end(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        sep = "|"
      )
      is_unique <- !duplicated(base_key)
      
    } else if (match_by == "stop_codon_transcript") {
      # Unique by stop position and transcript
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::end(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        as.character(GenomicRanges::mcols(base_gr)$transcript_id),
        sep = "|"
      )
      is_unique <- !duplicated(base_key)
      
    } else if (match_by == "start_stop") {
      # Unique by start and stop positions
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::start(base_gr),
        GenomicRanges::end(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        sep = "|"
      )
      is_unique <- !duplicated(base_key)
      
    } else if (match_by == "start_stop_transcript") {
      # Unique by start, stop, and transcript
      base_key <- paste(
        as.character(GenomicRanges::seqnames(base_gr)),
        GenomicRanges::start(base_gr),
        GenomicRanges::end(base_gr),
        as.character(GenomicRanges::strand(base_gr)),
        as.character(GenomicRanges::mcols(base_gr)$transcript_id),
        sep = "|"
      )
      is_unique <- !duplicated(base_key)
      
    } else if (match_by == "overlap") {
      # For overlap, need to group overlapping ORFs
      # Use reduce-like approach: iteratively merge overlapping ORFs
      hits <- GenomicRanges::findOverlaps(base_gr, base_gr)
      query_idx <- S4Vectors::subjectHits(hits)
      base_idx <- S4Vectors::queryHits(hits)
      
      overlap_width <- IRanges::width(GenomicRanges::pintersect(
        base_gr[base_idx], base_gr[query_idx]
      ))
      base_width <- IRanges::width(base_gr[base_idx])
      query_width <- IRanges::width(base_gr[query_idx])
      
      overlap_frac <- pmin(overlap_width / base_width,
                          overlap_width / query_width)
      
      # Build equivalence groups
      good_pairs <- hits[overlap_frac >= overlap_threshold]
      
      # Use first occurrence of each group
      groups <- rep(seq_along(base_gr), each = 1)
      for (i in seq_len(length(good_pairs))) {
        q <- S4Vectors::queryHits(good_pairs)[i]
        s <- S4Vectors::subjectHits(good_pairs)[i]
        groups[s] <- min(groups[q], groups[s])
      }
      is_unique <- !duplicated(groups)
      
    } else {  # metadata
      # Unique by metadata keys
      base_key <- do.call(paste, c(
        lapply(metadata_keys, function(k) {
          as.character(GenomicRanges::mcols(base_gr)[[k]])
        }),
        sep = "|"
      ))
      is_unique <- !duplicated(base_key)
    }
    
    base_gr <- base_gr[is_unique]
    presence_matrix <- presence_matrix[is_unique, , drop = FALSE]
  }
  
  # --- 5. Add metadata columns ------------------------------------------------
  # Retain only transcript_id and orf_type from the source data; all other
  # per-dataset metadata is dropped to avoid implying it represents the combined call.
  standard_cols <- c("transcript_id", "orf_type")
  if (keep_orf_ids) standard_cols <- c(standard_cols, "orf_id")
  keep_cols <- intersect(standard_cols, colnames(GenomicRanges::mcols(base_gr)))
  GenomicRanges::mcols(base_gr) <- GenomicRanges::mcols(base_gr)[, keep_cols, drop = FALSE]

  for (i in seq_len(n_datasets)) {
    col_name <- colnames(presence_matrix)[i]
    GenomicRanges::mcols(base_gr)[[col_name]] <- presence_matrix[, i]
  }

  # Add n_datasets count
  GenomicRanges::mcols(base_gr)$n_datasets <- as.integer(rowSums(presence_matrix))

  # Drop ORF IDs if requested
  if (!keep_orf_ids) {
    names(base_gr) <- NULL
  }
  
  base_gr
}


# ==============================================================================
# Export functions
# ==============================================================================

#' Export ORFs to BED12 with Spliced Exon Coordinates
#'
#' Maps parsed ORF predictions onto their canonical transcript's exon structure
#' and writes a BED12 file using \code{ORFik::export.bed12}. Each ORF is
#' represented as a multi-block entry that reflects the spliced exonic
#' structure of its \code{transcript_id}.
#'
#' Mapping is fully vectorised using \code{ORFik::pmapToTranscriptF} (genomic
#' to transcript space) and \code{ORFik::pmapFromTranscriptF} (transcript to
#' spliced genomic space), avoiding any per-ORF loop.
#'
#' @param orfs A GRanges object from \code{\link{parse_orfs}} (unified output).
#'   Must have a \code{transcript_id} metadata column.
#' @param txdb A TxDb object (e.g. from
#'   \code{GenomicFeatures::makeTxDbFromGFF()}) or a named GRangesList of
#'   exon ranges grouped by transcript (e.g. from
#'   \code{GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)}).
#'   Used to resolve the exon structure of each transcript.
#' @param file Output file path (should end in \code{.bed}).
#'
#' @return Invisibly returns the GRangesList of spliced ORFs written to disk.
#'
#' @details
#' The function first checks, vectorised via \code{findOverlaps}, whether each
#' ORF's start and stop positions fall within the exons of its assigned
#' \code{transcript_id}. For any ORF where they do not, it searches all
#' transcripts in \code{txdb} for an isoform whose exons cover both endpoints
#' and substitutes that isoform. If no compatible isoform exists, the ORF is
#' kept but exported as flat unspliced genomic coordinates (a single BED12
#' block). No ORFs are discarded. After re-assignment, a vectorised
#' element-wise \code{GRangesList::intersect} clips each ORF to the exon
#' blocks of its (possibly updated) transcript, and the result is written via
#' \code{ORFik::export.bed12}.
#'
#' ORFs whose original \code{transcript_id} is absent from \code{txdb} are
#' dropped with a warning (they cannot be mapped at all).
#'
#' @export
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomicRanges GRanges GRangesList findOverlaps intersect mcols
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @examples
#' \dontrun{
#' library(GenomicFeatures)
#' txdb <- makeTxDbFromGFF("gencode.v35.annotation.gtf")
#' orfs <- parse_orfs("ribotish_pred.txt", source = "ribotish")
#' export_orfs_bed(orfs, txdb, "ribotish_orfs_spliced.bed")
#' }
#'
#' @seealso \code{\link{parse_orfs}}, \code{ORFik::export.bed12},
#'   \code{ORFik::pmapToTranscriptF}, \code{ORFik::pmapFromTranscriptF}
export_orfs_bed <- function(orfs, txdb, file) {
  if (!is(orfs, "GRanges")) {
    stop("`orfs` must be a GRanges object from parse_orfs().")
  }
  if (!"transcript_id" %in% colnames(GenomicRanges::mcols(orfs))) {
    stop("`orfs` must have a `transcript_id` metadata column.")
  }
  if (!is.character(file) || length(file) != 1) {
    stop("`file` must be a single character string.")
  }

  # --- Resolve exon structure -------------------------------------------------
  if (is(txdb, "TxDb")) {
    tx_exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  } else if (is(txdb, "GRangesList")) {
    tx_exons <- txdb
  } else {
    stop("`txdb` must be a TxDb object or a named GRangesList.")
  }

  # --- Filter to transcripts present in txdb ----------------------------------
  tx_ids <- as.character(GenomicRanges::mcols(orfs)$transcript_id)
  found  <- tx_ids %in% names(tx_exons)

  if (!all(found)) {
    warning(sum(!found), " ORF(s) dropped: transcript_id not found in txdb.")
    orfs   <- orfs[found]
    tx_ids <- tx_ids[found]
  }

  if (length(orfs) == 0) {
    stop("No ORFs remain after filtering for known transcripts.")
  }

  # --- Resolve ORF names ------------------------------------------------------
  orf_names <- names(orfs)
  if (is.null(orf_names) || any(is.na(orf_names)) || any(orf_names == "")) {
    id_col    <- GenomicRanges::mcols(orfs)$orf_id
    orf_names <- if (!is.null(id_col)) as.character(id_col) else paste0("orf_", seq_along(orfs))
  }

  # Strip ORF metadata so intersection only sees coordinates
  orfs_coords <- orfs
  GenomicRanges::mcols(orfs_coords) <- NULL

  # --- Check that ORF start & stop fall within their assigned transcript ------
  # Build single-nucleotide markers for the ORF start and stop positions
  orf_starts_gr <- GenomicRanges::GRanges(
    GenomicRanges::seqnames(orfs_coords),
    IRanges::IRanges(GenomicRanges::start(orfs_coords), width = 1L),
    strand = GenomicRanges::strand(orfs_coords)
  )
  orf_stops_gr <- GenomicRanges::GRanges(
    GenomicRanges::seqnames(orfs_coords),
    IRanges::IRanges(GenomicRanges::end(orfs_coords), width = 1L),
    strand = GenomicRanges::strand(orfs_coords)
  )

  # Flatten the parallel transcript exons; track which element each row belongs to
  tx_par_init <- tx_exons[tx_ids]
  tx_par_flat <- unlist(tx_par_init, use.names = FALSE)
  tx_par_idx  <- rep(seq_along(tx_par_init), lengths(tx_par_init))

  hs_par <- GenomicRanges::findOverlaps(orf_starts_gr, tx_par_flat, ignore.strand = FALSE)
  he_par <- GenomicRanges::findOverlaps(orf_stops_gr,  tx_par_flat, ignore.strand = FALSE)

  # A hit is "self-consistent" when the ORF index equals the transcript element index
  valid_start <- seq_along(orfs_coords) %in%
    S4Vectors::queryHits(hs_par)[S4Vectors::queryHits(hs_par) ==
                                   tx_par_idx[S4Vectors::subjectHits(hs_par)]]
  valid_stop  <- seq_along(orfs_coords) %in%
    S4Vectors::queryHits(he_par)[S4Vectors::queryHits(he_par) ==
                                   tx_par_idx[S4Vectors::subjectHits(he_par)]]

  needs_remap <- !(valid_start & valid_stop)

  # --- Re-assign transcript for ORFs that extend beyond their assigned tx -----
  if (any(needs_remap)) {
    message(sum(needs_remap), " ORF(s) extend beyond their assigned transcript; ",
            "searching for a compatible isoform...")

    # Flatten the full tx_exons catalogue once
    tx_all_flat  <- unlist(tx_exons, use.names = FALSE)
    tx_all_names <- rep(names(tx_exons), lengths(tx_exons))

    remap_idx    <- which(needs_remap)
    remap_starts <- orf_starts_gr[remap_idx]
    remap_stops  <- orf_stops_gr[remap_idx]

    # Find all transcripts whose exons cover each re-mapped ORF's start/stop
    hs <- GenomicRanges::findOverlaps(remap_starts, tx_all_flat, ignore.strand = FALSE)
    he <- GenomicRanges::findOverlaps(remap_stops,  tx_all_flat, ignore.strand = FALSE)

    start_df <- data.frame(
      orf_i = S4Vectors::queryHits(hs),
      tx_id = tx_all_names[S4Vectors::subjectHits(hs)],
      stringsAsFactors = FALSE
    )
    stop_df <- data.frame(
      orf_i = S4Vectors::queryHits(he),
      tx_id = tx_all_names[S4Vectors::subjectHits(he)],
      stringsAsFactors = FALSE
    )

    # Keep only transcripts that cover both endpoints; take the first per ORF
    both_df <- merge(start_df, stop_df, by = c("orf_i", "tx_id"))
    best_df <- both_df[!duplicated(both_df$orf_i), ]

    # Apply re-assignments
    tx_ids[remap_idx[best_df$orf_i]] <- best_df$tx_id

    # ORFs with no compatible isoform: flag for genomic (unspliced) fallback
    genomic_fallback <- remap_idx[!seq_len(length(remap_idx)) %in% best_df$orf_i]
    if (length(genomic_fallback) > 0) {
      message(length(genomic_fallback),
              " ORF(s) could not be matched to any compatible isoform; ",
              "exporting as flat genomic (unspliced) coordinates.")
    }
  } else {
    genomic_fallback <- integer(0)
  }

  # --- Build parallel transcript GRangesList (with any re-assignments) --------
  tx_parallel <- tx_exons[tx_ids]

  # Convert flat GRanges to GRangesList (one element per ORF) for vectorised
  # element-wise intersection with tx_parallel
  orfs_grl <- as(orfs_coords, "GRangesList")

  # Vectorised element-wise intersection: keeps only exonic blocks of each ORF
  spliced_list <- GenomicRanges::intersect(orfs_grl, tx_parallel,
                                           ignore.strand = FALSE)

  # Replace fallback entries with a single-block genomic range
  if (length(genomic_fallback) > 0) {
    for (i in genomic_fallback) {
      spliced_list[[i]] <- orfs_coords[i]
    }
  }

  names(spliced_list) <- orf_names

  # --- Write BED12 via ORFik --------------------------------------------------
  ORFik::export.bed12(spliced_list, file)
  message("Exported ", length(spliced_list), " ORFs to ", file)
  invisible(spliced_list)
}


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Create unified output (GRanges with core metadata)
#' @keywords internal
.make_unified_output <- function(gr, additional_cols = NULL) {
  # Keep only essential metadata columns (plus any user-requested additional_cols)
  mcols_to_keep <- c("orf_id", "transcript_id", "orf_type", additional_cols)
  existing_cols <- mcols_to_keep[mcols_to_keep %in% colnames(GenomicRanges::mcols(gr))]
  
  if (length(existing_cols) > 0 && ncol(GenomicRanges::mcols(gr)) > length(existing_cols)) {
    GenomicRanges::mcols(gr) <- GenomicRanges::mcols(gr)[, existing_cols, drop = FALSE]
  }
  
  # Add orf_width (width of each range)
  GenomicRanges::mcols(gr)$orf_width <- IRanges::width(gr)
  
  gr
}

#' Create full output (list with unified GRangesList and metadata table)
#' @keywords internal
.make_full_output <- function(gr, raw_data) {
  # IMPORTANT: Extract metadata BEFORE calling .make_unified_output()
  # because .make_unified_output() strips metadata down to core columns
  
  # Extract all metadata into a data.frame
  # Prefer gr (processed data with mapped column names) over raw_data
  # except when raw_data is GRanges (e.g., ORFquant ORFs_tx with extra metadata)
  if (is(gr, "GRanges") || is(gr, "GRangesList")) {
    # Extract from gr first
    if (is(gr, "GRangesList")) {
      gr_unlisted <- unlist(gr, use.names = TRUE)
      metadata <- .mcols_to_dataframe(gr_unlisted)
    } else {
      metadata <- .mcols_to_dataframe(gr)
    }
  } else if (is(raw_data, "GRanges") || is(raw_data, "GRangesList")) {
    # If gr is not GRanges but raw_data is (e.g., ORFquant case),
    # extract from raw_data
    if (is(raw_data, "GRangesList")) {
      raw_gr <- unlist(raw_data, use.names = TRUE)
    } else {
      raw_gr <- raw_data
    }
    metadata <- .mcols_to_dataframe(raw_gr)
  } else if (is.data.frame(raw_data)) {
    # Last resort: use raw data.frame if gr is not GRanges
    metadata <- raw_data
  } else {
    # No metadata available
    metadata <- data.frame()
  }
  
  # Get unified output (grouped GRangesList)
  orfs <- .make_unified_output(gr, NULL)
  
  list(orfs = orfs, metadata = metadata)
}

#' Convert GRanges mcols to a simple data.frame
#' Handles complex column types like GRanges, AAStringSet, CompressedLists
#' @keywords internal
.mcols_to_dataframe <- function(gr) {
  mcols_df <- GenomicRanges::mcols(gr)
  
  # Add coordinate information first
  coord_df <- data.frame(
    seqnames = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    strand = as.character(GenomicRanges::strand(gr)),
    width = IRanges::width(gr),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  if (ncol(mcols_df) == 0) {
    return(coord_df)
  }
  
  # Fast path: Try direct conversion (works for simple types)
  metadata <- tryCatch({
    df <- as.data.frame(mcols_df, stringsAsFactors = FALSE)
    rownames(df) <- NULL
    df
  }, error = function(e) {
    # Slow path: Handle complex types one by one
    col_list <- list()
    for (i in seq_len(ncol(mcols_df))) {
      col_name <- colnames(mcols_df)[i]
      col_data <- mcols_df[[i]]
      
      col_list[[col_name]] <- tryCatch({
        if (is(col_data, "CompressedList") || is(col_data, "List")) {
          I(as.list(col_data))
        } else if (is(col_data, "AAStringSet")) {
          as.character(col_data)
        } else if (is(col_data, "GRanges") || is(col_data, "GRangesList")) {
          as.character(col_data)
        } else {
          col_data
        }
      }, error = function(e2) {
        tryCatch(as.character(col_data), 
                error = function(e3) I(as.list(col_data)))
      })
    }
    data.frame(col_list, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
  })
  
  # Combine with coordinates
  cbind(metadata, coord_df)
}

#' Read a file according to its spec
#' @keywords internal
.read_file <- function(file, spec) {

  # Custom reader takes priority
  if (!is.null(spec$read_fn)) {
    raw <- spec$read_fn(file)
    if (!is.data.frame(raw) && !is(raw, "GRanges") && !is(raw, "GRangesList")) {
      stop("read_fn must return a data.frame, GRanges, or GRangesList. ",
           "Got: ", class(raw)[1])
    }
    return(raw)
  }

  switch(spec$file_format,
    tsv = read.delim(file, stringsAsFactors = FALSE),
    csv = read.csv(file, stringsAsFactors = FALSE),
    bed = rtracklayer::import(file, format = "BED"),
    rds = readRDS(file),
    stop("Unknown file_format: ", spec$file_format)
  )
}


#' Resolve the first matching alias from a column_map entry
#'
#' @param canonical Canonical field name (e.g. "chrom").
#' @param column_map The spec's column_map list.
#' @param available Character vector of available column names in the data.
#' @param required Logical. If TRUE, error when no match found.
#' @return The matching column name, or NA if optional and not found.
#' @keywords internal
.resolve_column <- function(canonical, column_map, available,
                            required = FALSE) {
  if (!canonical %in% names(column_map)) {
    if (required) {
      stop("column_map has no entry for required field '", canonical, "'.")
    }
    return(NA_character_)
  }

  aliases <- column_map[[canonical]]
  hit <- intersect(aliases, available)

  if (length(hit) == 0) {
    if (required) {
      stop("Cannot find column for '", canonical,
           "'. Tried: ", paste(aliases, collapse = ", "),
           ". Available: ", paste(available, collapse = ", "))
    }
    return(NA_character_)
  }

  hit[1]
}


#' Build a GRanges from raw data and a spec
#' @keywords internal
.build_granges <- function(raw, spec, additional_cols = NULL) {

  # --- Already a GRanges (e.g. from BED import or RDS) ------------------------
  if (is(raw, "GRanges") || is(raw, "GRangesList")) {

    gr <- raw

    # For BED imports: set names from the 'name' column if available
    if (is.null(names(gr)) || all(names(gr) == "")) {
      if ("name" %in% colnames(GenomicRanges::mcols(gr))) {
        names(gr) <- GenomicRanges::mcols(gr)$name
      }
    }

    # Map canonical metadata if column_map has entries
    cols <- colnames(GenomicRanges::mcols(gr))
    cmap <- spec$column_map

    if ("orf_id" %in% names(cmap)) {
      col <- .resolve_column("orf_id", cmap, cols)
      if (!is.na(col)) {
        GenomicRanges::mcols(gr)$orf_id <- GenomicRanges::mcols(gr)[[col]]
        names(gr) <- GenomicRanges::mcols(gr)$orf_id
      }
    }
    if (is.null(GenomicRanges::mcols(gr)$orf_id)) {
      GenomicRanges::mcols(gr)$orf_id <- names(gr)
    }

    for (field in c("orf_type", "gene_id", "gene_name", "transcript_id")) {
      col <- .resolve_column(field, cmap, cols)
      if (!is.na(col) && !(field %in% colnames(GenomicRanges::mcols(gr)))) {
        GenomicRanges::mcols(gr)[[field]] <- GenomicRanges::mcols(gr)[[col]]
      }
    }

    # User-requested additional columns (for GRanges/GRangesList path)
    if (!is.null(additional_cols)) {
      for (col in additional_cols) {
        if (col %in% cols && !(col %in% colnames(GenomicRanges::mcols(gr)))) {
          GenomicRanges::mcols(gr)[[col]] <- GenomicRanges::mcols(gr)[[col]]
        } else if (!col %in% cols) {
          warning("Column '", col, "' not found in ", spec$name, " output.")
        }
      }
    }

    return(gr)
  }

  # --- data.frame path --------------------------------------------------------
  df <- raw
  cmap <- spec$column_map
  cols <- colnames(df)

  # Resolve required coordinate columns
  chr_col    <- .resolve_column("chrom",  cmap, cols, required = TRUE)
  start_col  <- .resolve_column("start",  cmap, cols, required = TRUE)
  end_col    <- .resolve_column("end",    cmap, cols, required = TRUE)
  strand_col <- .resolve_column("strand", cmap, cols, required = TRUE)
  id_col     <- .resolve_column("orf_id", cmap, cols, required = TRUE)

  start_pos <- as.integer(df[[start_col]])
  end_pos   <- as.integer(df[[end_col]])

  # Some callers (e.g. RiboCode) report start/stop in transcription direction,

  # so for minus-strand ORFs the genomic start > genomic stop.
  # GRanges requires start <= end, so swap where necessary.
  swap <- start_pos > end_pos
  if (any(swap)) {
    tmp <- start_pos[swap]
    start_pos[swap] <- end_pos[swap]
    end_pos[swap]   <- tmp
  }

  # Adjust 0-based starts (BED convention) -> 1-based (GRanges)
  if (spec$coord_system == "0-based") {
    start_pos <- start_pos + 1L
  }

  gr <- GenomicRanges::GRanges(
    seqnames = df[[chr_col]],
    ranges   = IRanges::IRanges(start = start_pos, end = end_pos),
    strand   = df[[strand_col]]
  )

  # Set names / orf_id
  orf_ids <- as.character(df[[id_col]])
  if (any(duplicated(orf_ids))) {
    orf_ids <- make.unique(orf_ids, sep = "_")
  }
  names(gr) <- orf_ids
  GenomicRanges::mcols(gr)$orf_id <- orf_ids

  # Map optional canonical metadata
  for (field in c("orf_type", "gene_id", "gene_name", "transcript_id")) {
    col <- .resolve_column(field, cmap, cols)
    if (!is.na(col)) {
      GenomicRanges::mcols(gr)[[field]] <- df[[col]]
    }
  }

  # Spec-declared extra columns
  for (col in spec$extra_cols) {
    if (col %in% cols && !(col %in% colnames(GenomicRanges::mcols(gr)))) {
      GenomicRanges::mcols(gr)[[col]] <- df[[col]]
    }
  }

  # User-requested additional columns
  if (!is.null(additional_cols)) {
    for (col in additional_cols) {
      if (col %in% cols && !(col %in% colnames(GenomicRanges::mcols(gr)))) {
        GenomicRanges::mcols(gr)[[col]] <- df[[col]]
      } else if (!col %in% cols) {
        warning("Column '", col, "' not found in ", spec$name, " output.")
      }
    }
  }

  gr
}


# ==============================================================================
# RiboCode-specific helpers
# ==============================================================================

#' Extend RiboCode ORF coordinates to include stop codon (simple genomic extension)
#'
#' RiboCode reports genomic coordinates excluding the stop codon.
#' This function extends each ORF by 3 nucleotides in genomic space
#' to include the stop codon (fast, but doesn't account for splicing).
#'
#' @param gr GRanges of ORFs
#' @param raw_df Raw RiboCode data.frame (not currently used)
#' @param txdb TxDb object (not used in simple version, kept for compatibility)
#' @return GRanges with extended coordinates
#' @keywords internal
.extend_ribocode_stop_codon <- function(gr, raw_df, txdb) {
  
  # Simple genomic extension:
  # Positive strand: decrease end by 3
  # Negative strand: increase start by 3
  is_plus <- as.character(GenomicRanges::strand(gr)) == "+"
  
  GenomicRanges::end(gr)[is_plus] <- GenomicRanges::end(gr)[is_plus] - 3L
  GenomicRanges::start(gr)[!is_plus] <- GenomicRanges::start(gr)[!is_plus] + 3L
  
  return(gr)
}


# ==============================================================================
# Validation and QC Functions
# ==============================================================================

#' Check Parsed ORF Caller Output
#'
#' Validates and summarizes the output from \code{\link{parse_orfs}}. Checks
#' for common issues like invalid coordinates, missing metadata, NA values,
#' and provides summary statistics.
#'
#' @param gr A GRanges object returned by \code{\link{parse_orfs}}.
#' @param verbose Logical. If TRUE (default), print detailed validation report.
#'   If FALSE, only return results invisibly.
#' @param check_width Logical. If TRUE (default), check that all ORF widths
#'   are multiples of 3 (expected for complete ORFs).
#' @param check_seqlevels Logical. If TRUE (default), check for non-standard
#'   chromosome names that might indicate genome version mismatches.
#'
#' @return Invisibly returns a list with validation results:
#'   \describe{
#'     \item{is_valid}{Logical. TRUE if all checks passed.}
#'     \item{n_orfs}{Total number of ORFs.}
#'     \item{issues}{Character vector of identified issues (empty if none).}
#'     \item{warnings}{Character vector of warnings (empty if none).}
#'     \item{summary}{Named list with summary statistics.}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' # Parse ORFs from RiboTISH
#' orfs <- parse_orfs("ribotish_output.txt", source = "ribotish")
#'
#' # Check the output
#' check_orf_calls(orfs)
#'
#' # Run silently and capture results
#' result <- check_orf_calls(orfs, verbose = FALSE)
#' if (!result$is_valid) {
#'   stop("ORF validation failed: ", paste(result$issues, collapse = "; "))
#' }
#' }
check_orf_calls <- function(gr, verbose = TRUE, check_width = TRUE, 
                           check_seqlevels = TRUE) {
  
  issues <- character(0)
  warnings <- character(0)
  
  # --- Basic Structure Checks ---
  if (!methods::is(gr, "GRanges")) {
    stop("Input must be a GRanges object. Got: ", class(gr)[1])
  }
  
  n_orfs <- length(gr)
  
  if (n_orfs == 0) {
    issues <- c(issues, "Empty GRanges (0 ORFs)")
  }
  
  # --- Coordinate Validation ---
  if (n_orfs > 0) {
    # Check for invalid ranges
    if (any(GenomicRanges::width(gr) < 1)) {
      n_invalid <- sum(GenomicRanges::width(gr) < 1)
      issues <- c(issues, sprintf("%d ORF(s) have width < 1", n_invalid))
    }
    
    # Check for unreasonably large ORFs (>50kb)
    if (any(GenomicRanges::width(gr) > 50000)) {
      n_large <- sum(GenomicRanges::width(gr) > 50000)
      warnings <- c(warnings, sprintf("%d ORF(s) are >50kb (unusually large)", n_large))
    }
    
    # Check for very small ORFs (<30nt = 10aa)
    if (any(GenomicRanges::width(gr) < 30)) {
      n_small <- sum(GenomicRanges::width(gr) < 30)
      warnings <- c(warnings, sprintf("%d ORF(s) are <30nt (very short)", n_small))
    }
    
    # Check if widths are multiples of 3
    if (check_width && any(GenomicRanges::width(gr) %% 3 != 0)) {
      n_not_mult3 <- sum(GenomicRanges::width(gr) %% 3 != 0)
      warnings <- c(warnings, sprintf("%d ORF(s) have length not divisible by 3", n_not_mult3))
    }
  }
  
  # --- Strand Validation ---
  if (n_orfs > 0) {
    strands <- as.character(GenomicRanges::strand(gr))
    n_unstranded <- sum(strands == "*")
    if (n_unstranded > 0) {
      warnings <- c(warnings, sprintf("%d ORF(s) have unspecified strand (*)", n_unstranded))
    }
    
    n_invalid_strand <- sum(!strands %in% c("+", "-", "*"))
    if (n_invalid_strand > 0) {
      issues <- c(issues, sprintf("%d ORF(s) have invalid strand", n_invalid_strand))
    }
  }
  
  # --- Seqlevel Checks ---
  if (n_orfs > 0 && check_seqlevels) {
    seqlevels <- GenomicRanges::seqlevels(gr)
    
    # Check for common issues
    has_chr_prefix <- any(grepl("^chr", seqlevels))
    has_no_chr_prefix <- any(!grepl("^chr", seqlevels) & !grepl("^scaffold|^contig", seqlevels))
    
    if (has_chr_prefix && has_no_chr_prefix) {
      warnings <- c(warnings, "Mixed chromosome naming (some with 'chr' prefix, some without)")
    }
    
    # Check for unusual names
    unusual <- seqlevels[grepl("^Un|_random|_alt|_fix", seqlevels)]
    if (length(unusual) > 0) {
      warnings <- c(warnings, sprintf("ORFs on %d unusual chromosome(s): %s", 
                                     length(unusual), paste(head(unusual, 3), collapse = ", ")))
    }
  }
  
  # --- Metadata Validation ---
  if (n_orfs > 0) {
    mcols_df <- as.data.frame(GenomicRanges::mcols(gr))
    
    if (ncol(mcols_df) == 0) {
      warnings <- c(warnings, "No metadata columns present")
    } else {
      # Check for NA values in key columns
      key_cols <- c("orf_id", "transcript_id", "gene_id", "orf_type")
      present_key_cols <- intersect(key_cols, colnames(mcols_df))
      
      for (col in present_key_cols) {
        n_na <- sum(is.na(mcols_df[[col]]))
        if (n_na > 0) {
          pct <- round(100 * n_na / n_orfs, 1)
          warnings <- c(warnings, sprintf("%d (%s%%) ORF(s) have NA in '%s'", n_na, pct, col))
        }
      }
      
      # Check for duplicate ORF IDs
      if ("orf_id" %in% colnames(mcols_df)) {
        orf_ids <- mcols_df$orf_id[!is.na(mcols_df$orf_id)]
        n_dup <- sum(duplicated(orf_ids))
        if (n_dup > 0) {
          issues <- c(issues, sprintf("%d duplicate ORF ID(s)", n_dup))
        }
      }
    }
  }
  
  # --- Summary Statistics ---
  summary_stats <- list(
    n_orfs = n_orfs,
    n_chromosomes = length(GenomicRanges::seqlevels(gr))
  )
  
  if (n_orfs > 0) {
    summary_stats$width_range <- range(GenomicRanges::width(gr))
    summary_stats$width_median <- median(GenomicRanges::width(gr))
    summary_stats$strand_counts <- table(as.character(GenomicRanges::strand(gr)))
    
    mcols_df <- as.data.frame(GenomicRanges::mcols(gr))
    if ("orf_type" %in% colnames(mcols_df)) {
      summary_stats$orf_type_counts <- table(mcols_df$orf_type, useNA = "ifany")
    }
    if ("transcript_id" %in% colnames(mcols_df)) {
      n_unique_tx <- length(unique(mcols_df$transcript_id[!is.na(mcols_df$transcript_id)]))
      summary_stats$n_unique_transcripts <- n_unique_tx
      summary_stats$orfs_per_transcript <- round(n_orfs / n_unique_tx, 2)
    }
  }
  
  # --- Determine Overall Status ---
  is_valid <- length(issues) == 0
  
  # --- Print Report ---
  if (verbose) {
    cat("\n")
    cat("========================================\n")
    cat("  ORF Caller Output Validation Report\n")
    cat("========================================\n\n")
    
    cat("Total ORFs:", n_orfs, "\n")
    if (n_orfs > 0) {
      cat("Chromosomes:", summary_stats$n_chromosomes, "\n")
      cat("Width range:", summary_stats$width_range[1], "-", summary_stats$width_range[2], "nt\n")
      cat("Median width:", summary_stats$width_median, "nt\n")
      
      cat("\nStrand distribution:\n")
      print(summary_stats$strand_counts)
      
      if (!is.null(summary_stats$orf_type_counts)) {
        cat("\nORF type distribution:\n")
        print(summary_stats$orf_type_counts)
      }
      
      if (!is.null(summary_stats$n_unique_transcripts)) {
        cat("\nUnique transcripts:", summary_stats$n_unique_transcripts, "\n")
        cat("ORFs per transcript:", summary_stats$orfs_per_transcript, "\n")
      }
    }
    
    cat("\n")
    if (is_valid && length(warnings) == 0) {
      cat("✓ All validation checks passed!\n")
    } else {
      if (length(issues) > 0) {
        cat("✗ ISSUES FOUND:\n")
        for (issue in issues) {
          cat("  •", issue, "\n")
        }
        cat("\n")
      }
      
      if (length(warnings) > 0) {
        cat("⚠ WARNINGS:\n")
        for (warn in warnings) {
          cat("  •", warn, "\n")
        }
        cat("\n")
      }
    }
    cat("========================================\n\n")
  }
  
  # --- Return Results ---
  result <- list(
    is_valid = is_valid,
    n_orfs = n_orfs,
    issues = issues,
    warnings = warnings,
    summary = summary_stats
  )
  
  invisible(result)
}


# ==============================================================================
# Built-in ORF Caller Specs
# ==============================================================================

#' Register Built-in ORF Caller Specs
#'
#' Called during package load to populate the registry with specs for
#' commonly used ORF callers. Users can override any built-in spec by
#' calling \code{register_orf_caller(..., overwrite = TRUE)}.
#'
#' @keywords internal
.register_builtin_specs <- function() {

  # ---------- Ribo-TISH -------------------------------------------------------
  register_orf_caller(orf_caller_spec(
    name        = "ribotish",
    description = "Ribo-TISH: Translation Initiation Site Hunter",
    file_format = "tsv",
    column_map  = list(
      chrom  = "chrom",
      start  = "start",
      end    = "end",
      strand = "strand",
      orf_id = "orf_id",
      orf_type      = "TisType",
      gene_name     = "Symbol",
      transcript_id = "Tid"
    ),
    coord_system = "1-based",
    url          = "https://github.com/zhpn1024/ribotish",
    extra_cols   = c("AALen", "StartCodon", "RiboPvalue"),
    read_fn = function(file) {
      df <- read.delim(file, stringsAsFactors = FALSE)
      # Ribo-TISH stores coordinates in GenomePos = "chr:start-end:strand"
      if ("GenomePos" %in% colnames(df)) {
        parts  <- strsplit(df$GenomePos, ":")
        df$chrom  <- vapply(parts, `[`, character(1), 1)
        coords    <- strsplit(vapply(parts, `[`, character(1), 2), "-")
        df$start  <- as.integer(vapply(coords, `[`, character(1), 1))
        df$end    <- as.integer(vapply(coords, `[`, character(1), 2))
        df$strand <- vapply(parts, `[`, character(1), 3)
      }
      # Build an orf_id if not present
      if (!"orf_id" %in% colnames(df)) {
        if ("Tid" %in% colnames(df)) {
          df$orf_id <- paste0("ribotish_", df$Tid, "_", df$start)
        } else {
          df$orf_id <- paste0("ribotish_orf_", seq_len(nrow(df)))
        }
      }
      df
    }
  ))

  # ---------- RiboCode --------------------------------------------------------
  register_orf_caller(orf_caller_spec(
    name        = "ribocode",
    description = "RiboCode: Detecting translated ORFs from Ribo-seq",
    file_format = "tsv",
    column_map  = list(
      chrom  = c("chrom", "Chrom"),
      start  = c("ORF_gstart", "orf_gstart"),
      end    = c("ORF_gstop", "orf_gstop"),
      strand = "strand",
      orf_id = c("ORF_ID", "orf_id"),
      orf_type      = c("ORF_type", "orf_type"),
      transcript_id = "transcript_id",
      gene_id       = "gene_id"
    ),
    coord_system = "1-based",
    url          = "https://github.com/xryanglab/RiboCode",
    extra_cols   = c("ORF_length", "pval"),
    post_process_fn = function(gr, raw_df, txdb = NULL, ribocode_extend_stop = TRUE) {
      # RiboCode excludes the stop codon.
      # We extend by 3 nt in genomic space to include it.
      
      # Skip extension if disabled
      if (!ribocode_extend_stop) {
        return(gr)
      }
      
      # Simple genomic extension (fast, no txdb required)
      .extend_ribocode_stop_codon(gr, raw_df, txdb)
    }
  ))

  # ---------- PRICE -----------------------------------------------------------
  register_orf_caller(orf_caller_spec(
    name        = "price",
    description = "PRICE: ORF Prediction from Ribosome Profiling",
    file_format = "bed",
    column_map  = list(
      orf_id = "name"
    ),
    coord_system = "0-based",
    url          = "https://github.com/erhard-lab/gedi"
  ))

  # ---------- Ribotricer ------------------------------------------------------
  register_orf_caller(orf_caller_spec(
    name        = "ribotricer",
    description = "Ribotricer: Detecting actively translating ORFs",
    file_format = "tsv",
    column_map  = list(
      chrom  = "chrom",
      start  = "start_codon",
      end    = "stop_codon",
      strand = "strand",
      orf_id = "ORF_ID",
      orf_type      = "ORF_type",
      transcript_id = "transcript_id"
    ),
    coord_system = "1-based",
    url          = "https://github.com/smithlabcode/ribotricer",
    extra_cols   = c("phase_score"),
    read_fn = function(file) {
      df <- read.delim(file, stringsAsFactors = FALSE)
      # If explicit coordinate columns are missing, extract from ORF_ID
      # Format: transcript__chrom__start__stop__strand__type
      if (!"start_codon" %in% colnames(df) && "ORF_ID" %in% colnames(df)) {
        id_parts <- strsplit(df$ORF_ID, "__")
        df$chrom       <- vapply(id_parts, `[`, character(1), 2)
        df$start_codon <- as.integer(vapply(id_parts, `[`, character(1), 3))
        df$stop_codon  <- as.integer(vapply(id_parts, `[`, character(1), 4))
        df$strand      <- vapply(id_parts, `[`, character(1), 5)
        if (!"ORF_type" %in% colnames(df)) {
          df$ORF_type <- vapply(id_parts, `[`, character(1), 6)
        }
      }
      df
    }
  ))

  # ---------- ORFquant --------------------------------------------------------

  # Helper: collapse a GRangesList (one element per ORF, each with >= 1 exon)
  # into a flat GRanges with one row per ORF using the full genomic span
  # (min start to max end). ORF-level metadata is taken from the first exon of
  # each element, since all exons share the same ORF-level attributes.
  # This ensures one row per ORF, consistent with single-range callers.
  .orfquant_grl_to_span <- function(grl) {
    span   <- unlist(range(grl), use.names = TRUE)
    n      <- S4Vectors::elementNROWS(grl)
    first  <- cumsum(n) - n + 1L
    flat   <- unlist(grl, use.names = FALSE)
    GenomicRanges::mcols(span) <- GenomicRanges::mcols(flat)[first, , drop = FALSE]
    GenomicRanges::mcols(span)$orf_id <- names(span)
    span
  }

  register_orf_caller(orf_caller_spec(
    name        = "orfquant",
    description = "ORFquant: Quantifying translation from Ribo-seq",
    file_format = "rds",
    column_map  = list(
      chrom  = c("Chromosome", "chromosome", "chrom"),
      start  = c("ORF_Start", "start"),
      end    = c("ORF_End", "end"),
      strand = c("Strand", "strand"),
      orf_id   = c("ORF_id_tr", "ORF_ID", "orf_id", "ORF_id"),
      orf_type = c("ORF_category_Tx", "ORF_category_Gen", "ORF_category", "ORF_type", "category"),
      gene_id  = "gene_id",
      gene_name = "gene_name",
      transcript_id = "transcript_id"
    ),
    coord_system = "1-based",
    url          = "https://github.com/lcalviell/ORFquant",
    extra_cols   = c("P_sites_raw", "P_sites_raw_uniq", "pval", 
                     "ORF_category_Tx", "ORF_category_Tx_compatible", "ORF_category_Gen"),
    post_process_fn = function(gr, raw_data, txdb = NULL, ribocode_extend_stop = TRUE) {
      # For ORFquant GRanges from ORFs_gen or ORFs_tx, ensure orf_type is set
      # Check both ORF_category_Tx and ORF_category_Gen
      gr_mcols <- GenomicRanges::mcols(gr)
      
      # Set orf_type if missing or all NA
      if (!"orf_type" %in% colnames(gr_mcols) || all(is.na(gr_mcols$orf_type))) {
        if ("ORF_category_Tx" %in% colnames(gr_mcols) && !all(is.na(gr_mcols$ORF_category_Tx))) {
          GenomicRanges::mcols(gr)$orf_type <- as.character(gr_mcols$ORF_category_Tx)
        } else if ("ORF_category_Gen" %in% colnames(gr_mcols) && !all(is.na(gr_mcols$ORF_category_Gen))) {
          GenomicRanges::mcols(gr)$orf_type <- as.character(gr_mcols$ORF_category_Gen)
        }
      }

      # ORFs_gen has one row per exon block; collapse to one row per ORF using
      # the full genomic span (min start to max end), carrying metadata from the
      # first exon of each ORF.
      orf_ids <- if ("orf_id" %in% colnames(GenomicRanges::mcols(gr))) {
        GenomicRanges::mcols(gr)$orf_id
      } else {
        names(gr)
      }
      if (!is.null(orf_ids) && any(duplicated(orf_ids))) {
        # S4Vectors::split uses PartitioningByEnd internally — faster than
        # base split() for large S4 GRanges objects.
        grl <- S4Vectors::split(gr, orf_ids)
        gr <- .orfquant_grl_to_span(grl)
      }

      gr
    },
    read_fn = function(file) {
      # ORFquant saves results as RData with ORFquant_results list.
      # Use ORFs_gen (genomic coordinates) and merge metadata from ORFs_tx if available.
      # Also supports plain RDS with GRanges or TSV.
      if (grepl("\\.rds$", file, ignore.case = TRUE)) {
        obj <- readRDS(file)
        if (is(obj, "GRanges") || is(obj, "GRangesList")) return(obj)
        if (is.data.frame(obj)) return(obj)
        if (is.list(obj) && "ORFs_gen" %in% names(obj)) return(obj$ORFs_gen)
        stop("ORFquant RDS must contain GRanges, GRangesList, data.frame, ",
             "or a list with ORFs_gen element.")
      }
      # Try loading as RData (e.g. final_ORFquant_results without extension)
      tryCatch({
        e <- new.env(parent = emptyenv())
        load(file, envir = e)
        obj_names <- ls(e)
        # Look for ORFquant_results list first
        if ("ORFquant_results" %in% obj_names) {
          res <- e$ORFquant_results
          
          # Strategy: Use ORFs_gen for genomic coordinates, merge metadata from ORFs_tx
          if ("ORFs_gen" %in% names(res) && "ORFs_tx" %in% names(res)) {
            # Both available - use genomic coords with transcript metadata
            gen <- res$ORFs_gen
            tx <- res$ORFs_tx
            
            # ORFs_gen may have no metadata columns but have names that correspond to ORF IDs
            # These names match ORF_id_tr in ORFs_tx
            gen_has_names <- !is.null(names(gen))
            gen_has_metadata <- ncol(GenomicRanges::mcols(gen)) > 0
            
            if (gen_has_names) {
              # Use names from ORFs_gen as orf_id and merge metadata from ORFs_tx
              GenomicRanges::mcols(gen)$orf_id <- names(gen)
              
              # Find the ORF ID column in ORFs_tx
              tx_mcols <- colnames(GenomicRanges::mcols(tx))
              orf_id_variants <- c("orf_id", "ORF_id", "ORF_ID", "ORF_id_tr", "ORF_id_gen")
              tx_orf_col <- intersect(orf_id_variants, tx_mcols)[1]
              
              if (!is.na(tx_orf_col)) {
                # Extract metadata columns from tx
                tx_metadata <- GenomicRanges::mcols(tx)
                metadata_cols <- setdiff(colnames(tx_metadata), 
                                         c("seqnames", "start", "end", "width", "strand"))
                
                # Match by orf_id
                gen_orf_ids <- names(gen)
                tx_orf_ids <- tx_metadata[[tx_orf_col]]
                match_idx <- match(gen_orf_ids, tx_orf_ids)
                
                # Bulk-assign all metadata columns at once via cbind on the DataFrame
                # to avoid repeated copy-on-modify from column-by-column assignment.
                cols_to_add <- intersect(metadata_cols, colnames(tx_metadata))
                if (length(cols_to_add) > 0) {
                  GenomicRanges::mcols(gen) <- cbind(
                    GenomicRanges::mcols(gen),
                    tx_metadata[match_idx, cols_to_add, drop = FALSE]
                  )
                }
                
                # Set orf_type from category columns
                gen_mcols_obj <- GenomicRanges::mcols(gen)
                if (!"orf_type" %in% colnames(gen_mcols_obj) || all(is.na(gen_mcols_obj$orf_type))) {
                  if ("ORF_category_Tx" %in% colnames(gen_mcols_obj) && !all(is.na(gen_mcols_obj$ORF_category_Tx))) {
                    GenomicRanges::mcols(gen)$orf_type <- as.character(gen_mcols_obj$ORF_category_Tx)
                  } else if ("ORF_category_Gen" %in% colnames(gen_mcols_obj) && !all(is.na(gen_mcols_obj$ORF_category_Gen))) {
                    GenomicRanges::mcols(gen)$orf_type <- as.character(gen_mcols_obj$ORF_category_Gen)
                  }
                }
                
                return(gen)
              }
            } else if (!gen_has_metadata) {
              # ORFs_gen has no names and no metadata - fall back to ORFs_tx
              warning("ORFquant ORFs_gen has no names or metadata. ",
                      "Using ORFs_tx (transcript coordinates) instead.")
              gen <- tx
              tx <- NULL  # Signal fallback
            }
            
            # Handle fallback case (tx is NULL)
            if (is.null(tx)) {
              # Handle standalone gen (which is actually ORFs_tx)
              if (is(gen, "GRangesList")) {
                gen <- .orfquant_grl_to_span(gen)
              } else {
                # Standardize orf_id column name
                gen_mcols_cols <- colnames(GenomicRanges::mcols(gen))
                orf_id_variants <- c("orf_id", "ORF_id", "ORF_ID", "ORF_id_tr", "ORF_id_gen")
                gen_orf_col <- intersect(orf_id_variants, gen_mcols_cols)[1]
                
                if (!is.na(gen_orf_col) && gen_orf_col != "orf_id") {
                  GenomicRanges::mcols(gen)$orf_id <- GenomicRanges::mcols(gen)[[gen_orf_col]]
                } else if (is.na(gen_orf_col)) {
                  if (!is.null(names(gen))) {
                    GenomicRanges::mcols(gen)$orf_id <- names(gen)
                  } else {
                    warning("No orf_id column found in ORFquant data, creating sequential IDs")
                    GenomicRanges::mcols(gen)$orf_id <- paste0("ORF_", seq_along(gen))
                  }
                }
              }
              
              # Set orf_type from category columns
              gen_mcols_obj <- GenomicRanges::mcols(gen)
              if (!"orf_type" %in% colnames(gen_mcols_obj) || all(is.na(gen_mcols_obj$orf_type))) {
                if ("ORF_category_Tx" %in% colnames(gen_mcols_obj) && !all(is.na(gen_mcols_obj$ORF_category_Tx))) {
                  GenomicRanges::mcols(gen)$orf_type <- as.character(gen_mcols_obj$ORF_category_Tx)
                } else if ("ORF_category_Gen" %in% colnames(gen_mcols_obj) && !all(is.na(gen_mcols_obj$ORF_category_Gen))) {
                  GenomicRanges::mcols(gen)$orf_type <- as.character(gen_mcols_obj$ORF_category_Gen)
                }
              }
              
              return(gen)
            }
            
            # If we reach here, we have both gen and tx with metadata columns to merge
            # (This handles the case where ORFs_gen has metadata in columns)
            
            # If GRangesList, collapse to one span per ORF
            if (is(gen, "GRangesList")) {
              gen <- .orfquant_grl_to_span(gen)
            }
            
            if (is(tx, "GRangesList")) {
              tx <- .orfquant_grl_to_span(tx)
            }
            
            # ORFs_gen has genomic coords but minimal metadata
            # ORFs_tx has transcript coords but rich metadata
            # Merge metadata from tx to gen based on orf_id (or variants)
            
            gen_mcols <- colnames(GenomicRanges::mcols(gen))
            tx_mcols <- colnames(GenomicRanges::mcols(tx))
            
            # Find orf_id column (try multiple possible names)
            orf_id_variants <- c("orf_id", "ORF_id", "ORF_ID", "ORF_id_tr", "ORF_id_gen")
            gen_orf_col <- intersect(orf_id_variants, gen_mcols)[1]
            tx_orf_col <- intersect(orf_id_variants, tx_mcols)[1]
            
            if (!is.na(gen_orf_col) && !is.na(tx_orf_col)) {
              # Extract metadata columns from tx (excluding coordinate-related)
              tx_metadata <- GenomicRanges::mcols(tx)
              metadata_cols <- setdiff(colnames(tx_metadata), 
                                       c("seqnames", "start", "end", "width", "strand"))
              
              # Match by orf_id
              gen_orf_ids <- GenomicRanges::mcols(gen)[[gen_orf_col]]
              tx_orf_ids <- tx_metadata[[tx_orf_col]]
              match_idx <- match(gen_orf_ids, tx_orf_ids)
              
              # Bulk-assign all metadata columns at once via cbind on the DataFrame
              # to avoid repeated copy-on-modify from column-by-column assignment.
              cols_to_add <- intersect(metadata_cols, colnames(tx_metadata))
              if (length(cols_to_add) > 0) {
                GenomicRanges::mcols(gen) <- cbind(
                  GenomicRanges::mcols(gen),
                  tx_metadata[match_idx, cols_to_add, drop = FALSE]
                )
              }
              
              # Ensure standard orf_id column exists (use the one from gen)
              if (gen_orf_col != "orf_id") {
                GenomicRanges::mcols(gen)$orf_id <- GenomicRanges::mcols(gen)[[gen_orf_col]]
              }
              
              # Ensure orf_type is set from category columns if not already present
              # Priority: ORF_category_Tx > ORF_category_Gen > existing orf_type
              gen_mcols_obj <- GenomicRanges::mcols(gen)
              if (!"orf_type" %in% colnames(gen_mcols_obj) || all(is.na(gen_mcols_obj$orf_type))) {
                if ("ORF_category_Tx" %in% colnames(gen_mcols_obj) && !all(is.na(gen_mcols_obj$ORF_category_Tx))) {
                  GenomicRanges::mcols(gen)$orf_type <- as.character(gen_mcols_obj$ORF_category_Tx)
                } else if ("ORF_category_Gen" %in% colnames(gen_mcols_obj) && !all(is.na(gen_mcols_obj$ORF_category_Gen))) {
                  GenomicRanges::mcols(gen)$orf_type <- as.character(gen_mcols_obj$ORF_category_Gen)
                }
              }
            } else {
              # No matching orf_id columns found - fallback to names or create sequential IDs
              if (!is.null(names(gen))) {
                # GRanges with names
                GenomicRanges::mcols(gen)$orf_id <- names(gen)
              } else {
                # Last resort: create sequential IDs
                warning("No orf_id column found in ORFquant data, creating sequential IDs")
                GenomicRanges::mcols(gen)$orf_id <- paste0("ORF_", seq_along(gen))
              }
            }
            return(gen)
          } else if ("ORFs_gen" %in% names(res)) {
            # Only genomic coords available - standardize orf_id column name
            gen <- res$ORFs_gen
            
            # If GRangesList, collapse to one span per ORF
            if (is(gen, "GRangesList")) {
              gen <- .orfquant_grl_to_span(gen)
            } else {
              # Handle regular GRanges
              gen_mcols <- colnames(GenomicRanges::mcols(gen))
              orf_id_variants <- c("orf_id", "ORF_id", "ORF_ID", "ORF_id_tr", "ORF_id_gen")
              gen_orf_col <- intersect(orf_id_variants, gen_mcols)[1]
              
              if (!is.na(gen_orf_col) && gen_orf_col != "orf_id") {
                GenomicRanges::mcols(gen)$orf_id <- GenomicRanges::mcols(gen)[[gen_orf_col]]
              } else if (is.na(gen_orf_col)) {
                # No standard orf_id column - try names or create IDs
                if (!is.null(names(gen))) {
                  GenomicRanges::mcols(gen)$orf_id <- names(gen)
                } else {
                  warning("No orf_id column found in ORFquant ORFs_gen, creating sequential IDs")
                  GenomicRanges::mcols(gen)$orf_id <- paste0("ORF_", seq_along(gen))
                }
              }
            }
            
            # Ensure orf_type is set from category columns
            gen_mcols_obj <- GenomicRanges::mcols(gen)
            if (!"orf_type" %in% colnames(gen_mcols_obj) || all(is.na(gen_mcols_obj$orf_type))) {
              if ("ORF_category_Gen" %in% colnames(gen_mcols_obj) && !all(is.na(gen_mcols_obj$ORF_category_Gen))) {
                GenomicRanges::mcols(gen)$orf_type <- as.character(gen_mcols_obj$ORF_category_Gen)
              } else if ("ORF_category_Tx" %in% colnames(gen_mcols_obj) && !all(is.na(gen_mcols_obj$ORF_category_Tx))) {
                GenomicRanges::mcols(gen)$orf_type <- as.character(gen_mcols_obj$ORF_category_Tx)
              }
            }
            
            return(gen)
          } else if ("ORFs_tx" %in% names(res)) {
            # Only transcript coords available - warn user
            warning("ORFquant file contains only ORFs_tx (transcript coordinates). ",
                    "Genomic coordinates (ORFs_gen) are preferred but not found.")
            
            tx <- res$ORFs_tx
            
            # If GRangesList, collapse to one span per ORF
            if (is(tx, "GRangesList")) {
              tx <- .orfquant_grl_to_span(tx)
            } else {
              # Handle regular GRanges
              tx_mcols <- colnames(GenomicRanges::mcols(tx))
              orf_id_variants <- c("orf_id", "ORF_id", "ORF_ID", "ORF_id_tr", "ORF_id_gen")
              tx_orf_col <- intersect(orf_id_variants, tx_mcols)[1]
              
              if (!is.na(tx_orf_col) && tx_orf_col != "orf_id") {
                GenomicRanges::mcols(tx)$orf_id <- GenomicRanges::mcols(tx)[[tx_orf_col]]
              } else if (is.na(tx_orf_col)) {
                # No standard orf_id column - try names or create IDs
                if (!is.null(names(tx))) {
                  GenomicRanges::mcols(tx)$orf_id <- names(tx)
                } else {
                  warning("No orf_id column found in ORFquant ORFs_tx, creating sequential IDs")
                  GenomicRanges::mcols(tx)$orf_id <- paste0("ORF_", seq_along(tx))
                }
              }
            }
            
            # Ensure orf_type is set from category columns
            tx_mcols_obj <- GenomicRanges::mcols(tx)
            if (!"orf_type" %in% colnames(tx_mcols_obj) || all(is.na(tx_mcols_obj$orf_type))) {
              if ("ORF_category_Tx" %in% colnames(tx_mcols_obj) && !all(is.na(tx_mcols_obj$ORF_category_Tx))) {
                GenomicRanges::mcols(tx)$orf_type <- as.character(tx_mcols_obj$ORF_category_Tx)
              } else if ("ORF_category_Gen" %in% colnames(tx_mcols_obj) && !all(is.na(tx_mcols_obj$ORF_category_Gen))) {
                GenomicRanges::mcols(tx)$orf_type <- as.character(tx_mcols_obj$ORF_category_Gen)
              }
            }
            
            return(tx)
          }
          stop("ORFquant_results list has no ORFs_gen or ORFs_tx element.")
        }
        # Otherwise look for any GRanges object
        for (nm in obj_names) {
          obj <- get(nm, envir = e)
          if (is(obj, "GRanges") || is(obj, "GRangesList")) return(obj)
        }
        stop("No GRanges found in ORFquant RData file.")
      }, error = function(err) {
        # Fall back to TSV
        tryCatch(
          read.delim(file, stringsAsFactors = FALSE),
          error = function(e2) stop("Cannot read ORFquant file: ", conditionMessage(err))
        )
      })
    }
  ))

  # ---------- Ribotie ---------------------------------------------------------
  register_orf_caller(orf_caller_spec(
    name        = "ribotie",
    description = "Ribotie: Ribo-seq ORF detection with TIS Transformer",
    file_format = "csv",
    column_map  = list(
      chrom         = "seqname",
      start         = "orf_start",
      end           = "orf_end",
      strand        = "strand",
      orf_id        = "ORF_id",
      orf_type      = "ORF_type",
      gene_id       = "gene_id",
      gene_name     = "gene_name",
      transcript_id = "transcript_id"
    ),
    coord_system = "1-based",
    url          = "https://github.com/TRISTAN-ORF/ribotie",
    extra_cols   = c("ribotie_score", "tis_transformer_score", "ribotie_rank",
                     "ORF_len", "start_codon", "reads_in_ORF",
                     "reads_in_frame_frac", "reads_coverage_frac"),
    read_fn = function(file) {
      df <- read.csv(file, stringsAsFactors = FALSE)
      # Ribotie provides TIS_coord (start) and LTS_coord (last coding nt)
      # and TTS_coord (first nt of stop codon). For GRanges we need the
      # ORF span: min(TIS, LTS) to max(TIS, LTS) on the correct strand.
      if (all(c("TIS_coord", "LTS_coord", "seqname", "strand") %in% colnames(df))) {
        tis <- as.integer(df$TIS_coord)
        lts <- as.integer(df$LTS_coord)
        df$orf_start <- pmin(tis, lts)
        df$orf_end   <- pmax(tis, lts)
      }
      df
    }
  ))

  # ---------- GENCODE Ribo-seq ORFs -------------------------------------------
  register_orf_caller(orf_caller_spec(
    name        = "gencode",
    description = "GENCODE: Consensus Ribo-seq ORF annotations",
    file_format = "bed",
    column_map  = list(
      orf_id = "name"
    ),
    coord_system = "0-based",
    url          = "https://www.gencodegenes.org/pages/riboseq_orfs/",
    post_process_fn = function(gr, raw_df, txdb = NULL, ribocode_extend_stop = TRUE) {
      GenomicRanges::mcols(gr)$orf_type <- "ribo-seq_orf"
      gr
    }
  ))
}


# -- Package load hook ---------------------------------------------------------

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  .register_builtin_specs()
}
