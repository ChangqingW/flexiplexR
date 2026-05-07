#' cutadapt wrapper
#' @description trim TSO adaptor with cutadapt
#' @param args arguments to be passed to cutadapt
#' @return Exit code of cutadapt
#'
#' @examples
#' cutadapt("-h")
#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskRun
#' @export
cutadapt <- function(args) {
  basiliskRun(env = flexiplexR_env, fun = function(x) {
    output <- base::system2("cutadapt", x, stdout = TRUE, stderr = "")
    return(output)
  }, x = args)
}
