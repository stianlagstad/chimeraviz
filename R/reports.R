#' Create a Fusion Report
#'
#' This function will create a html report with an overplot and a sortable,
#' searchable table with the fusion data.
#'
#' @param fusions A list of Fusion objects.
#' @param output_filename Output html-file filename.
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
#' fusions <- import_defuse(defuse833ke, "hg19", 3)
#' # Temporary file to store the report
#' output_filename <- tempfile(
#'   pattern = "fusionReport",
#'   fileext = ".html",
#'   tmpdir = tempdir())
#' # Create report
#' create_fusion_report(fusions, output_filename)
#'
#' @importFrom DT datatable
#' @importFrom plyr ldply
#'
#' @export
create_fusion_report <- function(
  fusions,
  output_filename,
  quiet = TRUE
) {

  # Check if we got a list of fusion objects
  if (class(fusions[[1]]) != "Fusion") {
    stop("fusions argument must be a list of Fusion objects")
  }

  # Get reference to the fusion list report file
  fusion_report_rmd <- system.file(
    "reports",
    "fusionListReport.Rmd",
    package = "chimeraviz")

  # Render the report
  rmarkdown::render(
    fusion_report_rmd,
    output_file = output_filename,
    params = list(fusions = fusions),
    quiet = quiet)
}
