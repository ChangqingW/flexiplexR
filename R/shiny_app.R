#' Launch the interactive barcode segment designer
#'
#' @description
#' Opens a Shiny web application for interactively building and reordering
#' barcode segment configurations. Segments can be dragged and dropped to
#' change their order. The app outputs the configuration as JSON (compatible
#' with FLAMES config files) and as R code (\code{barcode_segment()} calls).
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
#'
#' @return Invisibly returns \code{NULL}; the function is called for its
#'   side effect of launching the Shiny app.
#'
#' @examples
#' if (interactive()) {
#'   run_barcode_designer()
#' }
#'
#' @export
run_barcode_designer <- function(port = NULL, host = "127.0.0.1",
                                 launch.browser = FALSE) {
  app_dir <- system.file("shiny", "barcode_designer", package = "Rflexiplex")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Could not find barcode designer Shiny app in the Rflexiplex package.")
  }
  if (is.null(port)) {
    port <- httpuv::randomPort()
  }
  url <- sprintf("http://%s:%d", if (host == "0.0.0.0") "localhost" else host, port)
  message("Barcode designer running at: ", url)
  message("Press Ctrl+C to stop.")
  shiny::runApp(app_dir, port = port, host = host, launch.browser = launch.browser)
}
