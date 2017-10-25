#' Create a Fusion Report
#'
#' This function will create a html report with an overplot and a sortable,
#' searchable table with the fusion data.
#'
#' @param fusions A list of Fusion objects.
#' @param outputFilename Output html-file filename.
#' @param quiet Parameter passed to rmarkdown::render() to toggle its output.
#'
#' @return Creates a html report with an overplot and a sortable, searchable
#' table with the fusion data.
#'
#' @examples
#' # Load data
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 3)
#' # Temporary file to store the report
#' outputFilename <- tempfile(
#'   pattern = "fusionReport",
#'   fileext = ".html",
#'   tmpdir = tempdir())
#' # Create report
#' createFusionReport(fusions, outputFilename)
#'
#' @importFrom DT datatable
#' @importFrom plyr ldply
#'
#' @export
createFusionReport <- function(fusions, outputFilename, quiet = TRUE) {

  # Check if we got a list of fusion objects
  if (class(fusions[[1]]) != "Fusion") {
    stop("fusions argument must be a list of Fusion objects")
  }

  # Get reference to the fusion list report file
  fusionReportRmd <- system.file(
    "reports",
    "fusionListReport.Rmd",
    package="chimeraviz")

  # Render the report
  rmarkdown::render(
    fusionReportRmd,
    output_file = outputFilename,
    params = list(fusions = fusions),
    quiet = quiet)
}

createFusionReport2 <- function(fusions, outputFilename, ibam, edbFile, quiet = TRUE) {

  # Check if we got a list of fusion objects
  if (class(fusions[[1]]) != "Fusion") {
    stop("fusions argument must be a list of Fusion objects")
  }

  # Get reference to the fusion list report file
  fusionReportRmd <- system.file(
    "reports",
    "fusionListReport2.Rmd",
    package="chimeraviz")

  bam <- system.file(ibam)
  
  # Render the report
  rmarkdown::render(
    fusionReportRmd,
    output_file = outputFilename,
    params = list(fusions = fusions, bam = bam, edbFile = edbFile),
    quiet = quiet)
}

