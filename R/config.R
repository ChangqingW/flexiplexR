#' Read barcode parameters from a FLAMES-style JSON config file
#'
#' @description
#' Reads the \code{barcode_parameters} section from a FLAMES-style JSON
#' configuration file and returns it as a named list. The \code{segments}
#' element of the returned list is a data frame that can be passed directly
#' to \code{\link{find_barcode}} via its \code{segments} argument.
#'
#' @param config_file Path to a JSON configuration file containing a
#'   \code{barcode_parameters} top-level key.
#'
#' @return A named list with elements such as \code{segments} (data.frame),
#'   \code{barcode_groups}, \code{strand}, \code{TSO_seq},
#'   \code{TSO_prime}, \code{max_flank_editdistance}, etc.
#'
#' @examples
#' cfg_file <- system.file(
#'   "extdata", "config_sclr_nanopore_3end.json", package = "flexiplexR"
#' )
#' bp <- read_barcode_config(cfg_file)
#' str(bp$segments)
#'
#' @importFrom jsonlite fromJSON
#' @export
read_barcode_config <- function(config_file) {
  if (!file.exists(config_file)) {
    stop("Config file not found: ", config_file)
  }
  config <- jsonlite::fromJSON(config_file)
  bp <- config$barcode_parameters
  if (is.null(bp)) {
    stop("No 'barcode_parameters' section found in config: ", config_file)
  }
  bp
}
