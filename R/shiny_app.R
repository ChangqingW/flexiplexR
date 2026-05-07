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
#'   # With preview enabled using bundled fixtures
#'   bc <- tempfile()
#'   R.utils::gunzip(
#'     system.file("extdata", "bc_allow.tsv.gz", package = "Rflexiplex"),
#'     destname = bc, remove = FALSE
#'   )
#'   run_barcode_designer(
#'     fastq          = system.file("extdata", "fastq", "musc_rps24.fastq.gz",
#'                                  package = "Rflexiplex"),
#'     barcodes_files = bc
#'   )
#' }
#'
#' @export
run_barcode_designer <- function(port = NULL, host = "127.0.0.1",
                                 launch.browser = FALSE,
                                 fastq = NULL, barcodes_files = NULL) {
  app_dir <- system.file("shiny", "barcode_designer", package = "Rflexiplex")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Could not find barcode designer Shiny app in the Rflexiplex package.")
  }
  if (!is.null(fastq)) {
    stopifnot(length(fastq) == 1, is.character(fastq), file.exists(fastq))
  }
  if (!is.null(barcodes_files)) {
    stopifnot(is.character(barcodes_files), all(file.exists(barcodes_files)))
  }
  if (is.null(port)) {
    port <- httpuv::randomPort()
  }
  url <- sprintf("http://%s:%d", if (host == "0.0.0.0") "localhost" else host, port)
  message("Barcode designer running at: ", url)
  message("Press Ctrl+C to stop. Click 'Done' in the browser to return the configuration.")

  shiny::shinyOptions(
    rflexiplex_fastq          = fastq,
    rflexiplex_barcodes_files = barcodes_files
  )
  result <- shiny::runApp(app_dir, port = port, host = host, launch.browser = launch.browser)
  invisible(result)
}
