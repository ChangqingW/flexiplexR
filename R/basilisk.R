#' @importFrom basilisk BasiliskEnvironment
flexiplexR_env <- BasiliskEnvironment(
    envname = "flexiplexR_env", pkgname = "flexiplexR",
    packages = c(
        "python==3.11.9",
        "cutadapt==4.9"
    ),
    channels = c("conda-forge", "bioconda")
)
