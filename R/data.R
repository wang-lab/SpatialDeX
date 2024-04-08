#' This is an example dataset containing raw count data.
#'
#' @name example_raw_count
#' @title example_raw_count
#' @usage data(example_raw_count)
#' @format a data frame with 1000 rows (genes) and 50 variables (spots)
#' @keywords simulated spatial data
#' @source simulated raw count data for illustration purposes
#' @examples
#' data(example_raw_count)
"example_raw_count"
#'
#'
#' Annotation file for all known genes, sourced from copycat package
#'
#' @name full_anno
#' @title full_anno
#' @format a data frame with 56051 rows and 7 variables
#' @source this data is sourced from the `copycat` package
#'
full_anno <- load(system.file("extdata", "full_anno.RData", package = "SpatialDeX"))
#'
#'
#'
#' cyclegenes, sourced from copycat package
#'
#' @name cyclegenes
#' @title cyclegenes
#' @format a data frame with 1316 rows and 1 variables
#' @source this data is sourced from the `copycat` package
#'
cyclegenes <- load(system.file("extdata", "cyclegenes.RData", package = "SpatialDeX"))
#'
#' DNA.hg20, sourced from copycat package
#'
#' @name DNA.hg20
#' @title DNA.hg20
#' @format a data frame with 12205 rows and 3 variables
#' @source this data is sourced from the `copycat` package
#'
DNA.hg20 <- load(system.file("extdata", "DNA.hg20.RData", package = "SpatialDeX"))


