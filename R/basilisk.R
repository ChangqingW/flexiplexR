#' @importFrom basilisk BasiliskEnvironment
rflexiplex_env <- BasiliskEnvironment(
    envname = "rflexiplex_env", pkgname = "Rflexiplex",
    packages = c(
        "python==3.11.9",
        "cutadapt==4.9"
    ),
    channels = c("conda-forge", "bioconda")
)
