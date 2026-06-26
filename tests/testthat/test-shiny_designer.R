# Tests for the starting-configuration helpers behind run_barcode_designer().
# These exercise the normalization that converts the various accepted input
# shapes (JSON config, FlexiplexSegment objects, plain lists, data frames) into
# the plain-list format the Shiny app consumes.

bundled_cfg <- system.file("extdata", "config_sclr_nanopore_3end.json",
                           package = "flexiplexR")

test_that("FLAMES config file loads into plain segments with no NA leakage", {
  skip_if(bundled_cfg == "", "bundled config not installed")
  init <- flexiplexR:::.designer_init(config = bundled_cfg)

  expect_length(init$segments, 4)
  expect_identical(vapply(init$segments, `[[`, "", "type"),
                   c("FIXED", "MATCHED", "RANDOM", "FIXED"))
  expect_identical(vapply(init$segments, `[[`, "", "name"),
                   c("primer", "CB", "UB", "polyT"))

  # FIXED/RANDOM segments must carry "" (not NA) so nzchar() behaves in the app
  fixed <- init$segments[[1]]
  expect_identical(fixed$bc_list_name, "")
  expect_identical(fixed$group, "")
  expect_false(is.na(fixed$bc_list_name))

  matched <- init$segments[[2]]
  expect_identical(matched$bc_list_name, "CB")

  # numeric fields are integers, never NA
  bs <- vapply(init$segments, `[[`, 0L, "buffer_size")
  med <- vapply(init$segments, `[[`, 0L, "max_edit_distance")
  expect_true(is.integer(bs) && !anyNA(bs))
  expect_true(is.integer(med) && !anyNA(med))

  # global params come through from the config
  expect_identical(init$params$strand, "-")
  expect_equal(init$params$TSO_prime, 5)
  expect_identical(init$params$TSO_seq, "AAGCAGTGGTATCAACGCAGAGTACATGGG")

  # no barcode groups in this config
  expect_null(init$barcode_groups)
})

test_that("designer-exported JSON round-trips", {
  bp <- list(barcode_parameters = list(
    max_flank_editdistance  = 8,
    strand                  = "+",
    TSO_seq                 = "AGATCGGAAGAGCACA",
    TSO_prime               = 3,
    cutadapt_minimum_length = 10,
    full_length_only        = FALSE,
    segments = list(
      list(type = "FIXED",   pattern = "CTACAC", name = "primer",
           bc_list_name = "", group = "", buffer_size = 2L, max_edit_distance = 2L),
      list(type = "MATCHED", pattern = "NNNNNNNN", name = "CB",
           bc_list_name = "CB", group = "", buffer_size = 5L, max_edit_distance = 2L)
    ),
    barcode_groups = list()
  ))
  f <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(bp, pretty = TRUE, auto_unbox = TRUE), f)

  init <- flexiplexR:::.designer_init(config = f)
  expect_length(init$segments, 2)
  expect_identical(init$segments[[2]]$bc_list_name, "CB")
  expect_identical(init$segments[[2]]$buffer_size, 5L)
  expect_identical(init$params$strand, "+")
  expect_equal(init$params$TSO_prime, 3)
  expect_null(init$barcode_groups)
})

test_that("FlexiplexSegment objects normalize to plain lists", {
  segs <- list(
    barcode_segment("FIXED",   "CTACAC", "primer"),
    barcode_segment("MATCHED", "NNNN", "CB", bc_list = "CB",
                    buffer_size = 5, max_edit_distance = 2),
    barcode_segment("RANDOM",  "NNNN", "UB")
  )
  init <- flexiplexR:::.designer_init(segments = segs)

  expect_length(init$segments, 3)
  expect_identical(init$segments[[1]]$bc_list_name, "")   # FIXED -> "" not NA
  expect_identical(init$segments[[3]]$group, "")          # RANDOM -> "" not NA
  expect_identical(init$segments[[2]]$bc_list_name, "CB")
  expect_identical(init$segments[[2]]$buffer_size, 5L)
  expect_null(init$params)                                # lists carry no params
})

test_that("data.frame from read_barcode_config normalizes", {
  skip_if(bundled_cfg == "", "bundled config not installed")
  df <- read_barcode_config(bundled_cfg)$segments
  expect_true(is.data.frame(df))

  init <- flexiplexR:::.designer_init(segments = df)
  expect_length(init$segments, 4)
  expect_identical(init$segments[[1]]$bc_list_name, "")
  expect_identical(init$segments[[2]]$bc_list_name, "CB")
})

test_that("barcode groups normalize and empty groups stay empty", {
  grps <- list(barcode_group("CB", "CB", max_edit_distance = 3))
  out <- flexiplexR:::.normalize_designer_groups(grps)
  expect_length(out, 1)
  expect_identical(out[[1]]$name, "CB")
  expect_identical(out[[1]]$bc_list_name, "CB")
  expect_identical(out[[1]]$max_edit_distance, 3L)

  expect_null(flexiplexR:::.normalize_designer_groups(list()))
  expect_null(flexiplexR:::.normalize_designer_groups(NULL))
})

test_that("config and segments are mutually exclusive", {
  skip_if(bundled_cfg == "", "bundled config not installed")
  expect_error(
    flexiplexR:::.designer_init(
      config   = bundled_cfg,
      segments = list(barcode_segment("FIXED", "ACGT", "primer"))
    ),
    "not both"
  )
})

test_that("no input falls through to NULL (app uses defaults)", {
  init <- flexiplexR:::.designer_init()
  expect_null(init$segments)
  expect_null(init$barcode_groups)
  expect_null(init$params)
})

test_that(".cutadapt_fulllength computes the full-length rate", {
  rep <- list(read_counts = list(input = 200L, output = 150L,
                                 read1_with_adapter = 120L,
                                 read2_with_adapter = NULL))
  fl <- flexiplexR:::.cutadapt_fulllength(rep)
  expect_identical(fl$input, 200L)
  expect_identical(fl$with_adapter, 120L)
  expect_identical(fl$output, 150L)
  expect_equal(fl$rate, 60)
})

test_that(".cutadapt_fulllength returns NULL when cutadapt did not run", {
  expect_null(flexiplexR:::.cutadapt_fulllength(NULL))
  expect_null(flexiplexR:::.cutadapt_fulllength(list()))          # no read_counts
  expect_null(flexiplexR:::.cutadapt_fulllength(list(read_counts = list())))
})

test_that(".cutadapt_fulllength is divide-by-zero safe and tolerates missing fields", {
  fl <- flexiplexR:::.cutadapt_fulllength(
    list(read_counts = list(input = 0L, read1_with_adapter = 0L))
  )
  expect_equal(fl$rate, 0)
  expect_true(is.na(fl$output))

  # missing read1_with_adapter -> treated as 0
  fl2 <- flexiplexR:::.cutadapt_fulllength(
    list(read_counts = list(input = 10L, output = 10L))
  )
  expect_identical(fl2$with_adapter, 0L)
  expect_equal(fl2$rate, 0)
})
