#' Launch the interactive barcode segment designer
#'
#' @description
#' Opens a Shiny web application for interactively building and reordering
#' barcode segment configurations. Segments can be dragged and dropped to
#' change their order. The app outputs the configuration as JSON (compatible
#' with FLAMES config files) and as R code (\code{barcode_segment()} calls).
#'
#' If \code{fastq} and \code{barcodes_files} are supplied, a Preview tab is
#' enabled that runs \code{\link{find_barcode}} on a small subsample of the
#' input and reports the demultiplex rate together with example demultiplexed
#' read headers (which carry the \code{CB:Z} / \code{UB:Z} tags) and rows from
#' the stats output. The preview is omitted when either argument is missing.
#'
#' On headless servers or SSH sessions, set \code{launch.browser = FALSE}
#' (the default when no display is detected) and open the printed URL in
#' your local browser, optionally via SSH port forwarding:
#' \preformatted{ssh -L <port>:localhost:<port> user@host}
#'
#' @param port Port to listen on. Defaults to a random available port.
#' @param host Host address to bind. Use \code{"0.0.0.0"} to allow remote
#'   connections (e.g. over SSH tunnel). Defaults to \code{"127.0.0.1"}.
#' @param launch.browser Whether to open a browser automatically. Defaults to
#'   \code{FALSE}; open the printed URL in your local browser instead, using
#'   SSH port forwarding if needed:
#'   \preformatted{ssh -L <port>:localhost:<port> user@host}
#' @param fastq Optional path to a FASTQ (or gzipped FASTQ) file used to
#'   populate the Preview tab. When \code{NULL} (default) the Preview tab is
#'   hidden.
#' @param barcodes_files Optional barcode allow-list path(s) for the Preview
#'   tab. Either a single unnamed character string (used for all
#'   \code{MATCHED} segments) or a named character vector keyed by
#'   \code{bc_list_name}. Same shape as the \code{barcodes_files} argument of
#'   \code{\link{find_barcode}}.
#' @param config Optional starting configuration. Either a path to a JSON file
#'   or an already-parsed list. JSON may be in the format this designer itself
#'   downloads \emph{or} a full FLAMES config (both wrap the relevant settings
#'   in a \code{barcode_parameters} key); it is read with
#'   \code{\link{read_barcode_config}}. A list may be a parsed
#'   \code{barcode_parameters} list or the \code{list(segments, barcode_groups)}
#'   value returned by this function's Done button. When supplied, the global
#'   parameters (strand, TSO, edit distances, etc.) present in the config are
#'   also pre-filled. Mutually exclusive with \code{segments}/\code{barcode_groups}.
#' @param segments Optional starting segments, used when \code{config} is not
#'   given. A list of \code{\link{barcode_segment}} objects, a list of plain
#'   lists (as returned by the Done button), or a data frame (as returned by
#'   \code{\link{read_barcode_config}}'s \code{segments} element). Global
#'   parameters keep their defaults on this path.
#' @param barcode_groups Optional starting barcode groups, paired with
#'   \code{segments}. A list of \code{\link{barcode_group}} objects, a list of
#'   plain lists, or a data frame. Defaults to no groups.
#'
#' @return Invisibly returns the value passed to \code{shiny::stopApp()} by
#'   the Done button: a list with \code{segments} and \code{barcode_groups}.
#'   Returns \code{NULL} if the app is closed without clicking Done.
#'
#' @examples
#' if (interactive()) {
#'   # Designer only
#'   run_barcode_designer()
#'
#'   # Start from an existing config (FLAMES or designer-exported JSON)
#'   run_barcode_designer(
#'     config = system.file("extdata", "config_sclr_nanopore_3end.json",
#'                          package = "flexiplexR")
#'   )
#'
#'   # Start from segment objects directly
#'   run_barcode_designer(
#'     segments = list(
#'       barcode_segment("FIXED",   "CTACACGACGCTCTTCCGATCT", "primer"),
#'       barcode_segment("MATCHED", "NNNNNNNNNNNNNNNN", "CB", bc_list = "CB"),
#'       barcode_segment("RANDOM",  "NNNNNNNNNNNN", "UB")
#'     )
#'   )
#'
#'   # With preview enabled using bundled fixtures
#'   bc <- tempfile()
#'   R.utils::gunzip(
#'     system.file("extdata", "bc_allow.tsv.gz", package = "flexiplexR"),
#'     destname = bc, remove = FALSE
#'   )
#'   run_barcode_designer(
#'     fastq          = system.file("extdata", "fastq", "musc_rps24.fastq.gz",
#'                                  package = "flexiplexR"),
#'     barcodes_files = bc
#'   )
#' }
#'
#' @export
run_barcode_designer <- function(port = NULL, host = "127.0.0.1",
                                 launch.browser = FALSE,
                                 fastq = NULL, barcodes_files = NULL,
                                 config = NULL,
                                 segments = NULL, barcode_groups = NULL) {
  app_dir <- system.file("shiny", "barcode_designer", package = "flexiplexR")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Could not find barcode designer Shiny app in the flexiplexR package.")
  }
  if (!is.null(fastq)) {
    stopifnot(length(fastq) == 1, is.character(fastq), file.exists(fastq))
  }
  if (!is.null(barcodes_files)) {
    stopifnot(is.character(barcodes_files), all(file.exists(barcodes_files)))
  }

  init <- .designer_init(config, segments, barcode_groups)

  if (is.null(port)) {
    port <- httpuv::randomPort()
  }
  url <- sprintf("http://%s:%d", if (host == "0.0.0.0") "localhost" else host, port)
  message("Barcode designer running at: ", url)
  message("Press Ctrl+C to stop. Click 'Done' in the browser to return the configuration.")

  shiny::shinyOptions(
    flexiplexR_fastq          = fastq,
    flexiplexR_barcodes_files = barcodes_files,
    flexiplexR_init_segments  = init$segments,
    flexiplexR_init_groups    = init$barcode_groups,
    flexiplexR_init_params    = init$params
  )
  result <- shiny::runApp(app_dir, port = port, host = host, launch.browser = launch.browser)
  invisible(result)
}

# Resolve the various accepted starting-configuration shapes into the plain-list
# format the Shiny app consumes. Returns a list with `segments` (list of plain
# lists or NULL), `barcode_groups` (list of plain lists or NULL) and `params`
# (named list of global parameters or NULL).
.designer_init <- function(config = NULL, segments = NULL, barcode_groups = NULL) {
  if (!is.null(config) && (!is.null(segments) || !is.null(barcode_groups))) {
    stop("Provide either `config` or `segments`/`barcode_groups`, not both.")
  }

  if (is.null(config)) {
    return(list(
      segments       = .normalize_designer_segments(segments),
      barcode_groups = .normalize_designer_groups(barcode_groups),
      params         = NULL
    ))
  }

  if (is.character(config)) {
    stopifnot(length(config) == 1)
    bp <- read_barcode_config(config)
  } else if (is.list(config)) {
    bp <- config
  } else {
    stop("`config` must be a JSON file path or a list.")
  }

  param_names <- c("max_flank_editdistance", "strand", "TSO_seq", "TSO_prime",
                   "cutadapt_minimum_length", "full_length_only")
  params <- bp[intersect(param_names, names(bp))]
  if (length(params) == 0) params <- NULL

  list(
    segments       = .normalize_designer_segments(bp$segments),
    barcode_groups = .normalize_designer_groups(bp$barcode_groups),
    params         = params
  )
}

# Read a named field from a segment/group regardless of whether it is an S4
# FlexiplexSegment/FlexiplexGroup object or a plain list / data.frame row.
# NULL and NA collapse to `default`.
.designer_field <- function(x, field, default) {
  v <- if (isS4(x)) {
    if (field %in% methods::slotNames(x)) methods::slot(x, field) else default
  } else {
    x[[field]]
  }
  if (is.null(v) || (length(v) == 1 && is.na(v))) default else v
}

# Normalize a segments input (list of FlexiplexSegment objects, list of plain
# lists, or a data.frame) into a list of plain lists with all seven fields.
.normalize_designer_segments <- function(segs) {
  if (is.null(segs) || length(segs) == 0) return(NULL)
  if (is.data.frame(segs)) {
    segs <- lapply(seq_len(nrow(segs)), function(i) as.list(segs[i, , drop = FALSE]))
  }
  lapply(segs, function(s) {
    type <- as.character(.designer_field(s, "type", "FIXED"))
    list(
      type              = type,
      pattern           = as.character(.designer_field(s, "pattern", "")),
      name              = as.character(.designer_field(s, "name", type)),
      bc_list_name      = as.character(.designer_field(s, "bc_list_name", "")),
      group             = as.character(.designer_field(s, "group", "")),
      buffer_size       = as.integer(.designer_field(s, "buffer_size", 2L)),
      max_edit_distance = as.integer(.designer_field(s, "max_edit_distance", 2L))
    )
  })
}

# Normalize a barcode_groups input (list of FlexiplexGroup objects, list of
# plain lists, or a data.frame) into a list of plain lists.
.normalize_designer_groups <- function(grps) {
  if (is.null(grps) || length(grps) == 0) return(NULL)
  if (is.data.frame(grps)) {
    grps <- lapply(seq_len(nrow(grps)), function(i) as.list(grps[i, , drop = FALSE]))
  }
  lapply(grps, function(g) {
    list(
      name              = as.character(.designer_field(g, "name", "")),
      bc_list_name      = as.character(.designer_field(g, "bc_list_name", "")),
      max_edit_distance = as.integer(.designer_field(g, "max_edit_distance", 2L))
    )
  })
}

# Summarize cutadapt's adapter-finding stats into a full-length-read estimate.
# `cutadapt` is the parsed cutadapt JSON report returned by find_barcode() as
# report$cutadapt, or NULL when no TSO / adapter sequence was supplied (then this
# returns NULL). cutadapt searches for the far-end adapter, so the fraction of
# reads where it is found estimates the proportion of *full-length* reads. The
# denominator `input` is the number of reads fed to cutadapt, i.e. the
# demultiplexed reads. (`read1_with_adapter` is the cutadapt 4.9 read_counts
# field — see R/basilisk.R for the pinned version; this is the single place that
# field name is referenced, so version drift is a one-line fix here.)
.cutadapt_fulllength <- function(cutadapt) {
  rc <- cutadapt$read_counts
  if (is.null(rc) || is.null(rc$input)) return(NULL)
  input        <- as.integer(rc$input)
  with_adapter <- rc$read1_with_adapter
  with_adapter <- if (is.null(with_adapter) || is.na(with_adapter)) 0L else as.integer(with_adapter)
  output       <- if (is.null(rc$output)) NA_integer_ else as.integer(rc$output)
  rate <- if (length(input) == 1 && !is.na(input) && input > 0)
            round(100 * with_adapter / input, 2) else 0
  list(input = input, with_adapter = with_adapter, output = output, rate = rate)
}
