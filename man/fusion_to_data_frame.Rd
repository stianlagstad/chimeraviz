% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{fusion_to_data_frame}
\alias{fusion_to_data_frame}
\title{Coerce Fusion object to data.frame}
\usage{
fusion_to_data_frame(fusion)
}
\arguments{
\item{fusion}{The Fusion object to coerce.}
}
\value{
A data.frame with the fusion object.
}
\description{
This function is used in create_fusion_report() to convert Fusion objects to a
data.frame-format.
}
\examples{
# Load data
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 1)
# Find the fusion object to create a data frame from
fusion <- get_fusion_by_id(fusions, 5267)
# Create the data frame
dfFusion <- fusion_to_data_frame(fusion)

}
\seealso{
create_fusion_report
}
