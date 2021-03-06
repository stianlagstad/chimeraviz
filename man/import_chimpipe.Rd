% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_chimpipe.R
\name{import_chimpipe}
\alias{import_chimpipe}
\title{Import results from a ChimPipe run into a list of Fusion objects.}
\usage{
import_chimpipe(filename, genome_version, limit)
}
\arguments{
\item{filename}{Filename for the ChimPipe results.}

\item{genome_version}{Which genome was used in mapping (hg19, hg38, etc.).}

\item{limit}{A limit on how many lines to read.}
}
\value{
A list of Fusion objects.
}
\description{
A function that imports the results from a ChimPipe run, typically from a
chimericJunctions_.txt file, into a list of Fusion objects.
}
\examples{
chimpipefile <- system.file(
  "extdata",
  "chimericJunctions_MCF-7.txt",
  package="chimeraviz")
fusions <- import_chimpipe(chimpipefile, "hg19", 3)
# This should import a list of 3 fusions described in Fusion objects.
}
