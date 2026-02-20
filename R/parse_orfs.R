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
#'   to retain. Only used when \code{output = "full"}.
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
                       additional_cols = NULL) {

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
    if (is.data.frame(raw)) {
      gr <- spec$post_process_fn(gr, raw)
    } else {
      gr <- spec$post_process_fn(gr, NULL)
    }
  }

  # --- 4.5. Strip metadata early for unified output (performance optimization) --
  # For unified output, we only need transcript_id, orf_type, orf_id
  # Stripping unnecessary columns (especially complex types like CompressedLists)
  # before genome style conversion significantly improves performance
  if (output == "unified") {
    cols_to_keep <- c("transcript_id", "orf_type", "orf_id")
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
    return(.make_unified_output(gr))
  } else {
    return(.make_full_output(gr, raw))
  }
}


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Create unified output (GRanges with core metadata)
#' @keywords internal
.make_unified_output <- function(gr) {
  # Keep only essential metadata columns
  mcols_to_keep <- c("orf_id", "transcript_id", "orf_type")
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
  orfs <- .make_unified_output(gr)
  
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
    extra_cols   = c("ORF_length", "pval")
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
    post_process_fn = function(gr, raw_data) {
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
                
                # Copy metadata columns directly from tx to gen
                for (col in metadata_cols) {
                  if (col %in% colnames(tx_metadata)) {
                    GenomicRanges::mcols(gen)[[col]] <- tx_metadata[[col]][match_idx]
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
                orf_ids_gen <- rep(names(gen), elementNROWS(gen))
                gen <- unlist(gen, use.names = FALSE)
                GenomicRanges::mcols(gen)$orf_id <- orf_ids_gen
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
            
            # If GRangesList, unlist and propagate names as orf_id
            if (is(gen, "GRangesList")) {
              orf_ids_gen <- rep(names(gen), elementNROWS(gen))
              gen <- unlist(gen, use.names = FALSE)
              if (is.null(GenomicRanges::mcols(gen)$orf_id)) {
                GenomicRanges::mcols(gen)$orf_id <- orf_ids_gen
              }
            }
            
            if (is(tx, "GRangesList")) {
              orf_ids_tx <- rep(names(tx), elementNROWS(tx))
              tx <- unlist(tx, use.names = FALSE)
              if (is.null(GenomicRanges::mcols(tx)$orf_id)) {
                GenomicRanges::mcols(tx)$orf_id <- orf_ids_tx
              }
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
              
              # Directly copy metadata columns from tx to gen (avoid data.frame conversion for complex types)
              for (col in metadata_cols) {
                if (col %in% colnames(tx_metadata)) {
                  # Copy the column directly, preserving complex types
                  GenomicRanges::mcols(gen)[[col]] <- tx_metadata[[col]][match_idx]
                }
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
            
            # If GRangesList, unlist and propagate names as orf_id
            if (is(gen, "GRangesList")) {
              orf_ids_gen <- rep(names(gen), elementNROWS(gen))
              gen <- unlist(gen, use.names = FALSE)
              GenomicRanges::mcols(gen)$orf_id <- orf_ids_gen
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
            
            # If GRangesList, unlist and propagate names as orf_id
            if (is(tx, "GRangesList")) {
              orf_ids_tx <- rep(names(tx), elementNROWS(tx))
              tx <- unlist(tx, use.names = FALSE)
              GenomicRanges::mcols(tx)$orf_id <- orf_ids_tx
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
    post_process_fn = function(gr, raw_df) {
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
