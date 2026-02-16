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
#' @param additional_cols Character vector of extra caller-specific columns
#'   to retain as GRanges metadata.
#'
#' @return A named GRanges object with canonical metadata columns
#'   (\code{orf_id}, and optionally \code{orf_type}, \code{gene_id},
#'   \code{gene_name}, \code{transcript_id}).
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
#' }
#'
#' @export
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @importFrom GenomicRanges GRanges mcols mcols<-
#' @importFrom IRanges IRanges
#'
#' @examples
#' # Parse GENCODE consensus ORFs (ships with package)
#' bed <- system.file("extdata", "Ribo-seq_ORFs.bed", package = "BilbORF")
#' if (bed != "") {
#'   orfs <- parse_orfs(bed, source = "gencode")
#'   print(length(orfs))
#' }
#'
#' \dontrun{
#' # Parse Ribo-TISH results
#' orfs <- parse_orfs("ribotish_pred.txt", source = "ribotish")
#'
#' # Parse RiboCode with Ensembl chromosomes
#' orfs <- parse_orfs("ribocode.txt",
#'                    source = "ribocode",
#'                    genome_style = "Ensembl")
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
                       additional_cols = NULL) {

  # --- Validate inputs --------------------------------------------------------
  if (!file.exists(file)) stop("File not found: ", file)
  if (!genome_style %in% c("UCSC", "Ensembl")) {
    stop("genome_style must be 'UCSC' or 'Ensembl'.")
  }

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

  gr
}


# ==============================================================================
# Internal helpers
# ==============================================================================

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
      orf_type = c("ORF_category_Tx", "ORF_category", "ORF_type", "category"),
      gene_id  = "gene_id",
      gene_name = "gene_name",
      transcript_id = "transcript_id"
    ),
    coord_system = "1-based",
    url          = "https://github.com/lcalviell/ORFquant",
    extra_cols   = c("P_sites_raw", "P_sites_raw_uniq", "pval"),
    read_fn = function(file) {
      # ORFquant saves results as RData with ORFquant_results list.
      # The key element is ORFs_gen (genomic GRanges) or ORFs_tx (tx-level).
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
          if ("ORFs_gen" %in% names(res)) return(res$ORFs_gen)
          if ("ORFs_tx" %in% names(res)) return(res$ORFs_tx)
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
