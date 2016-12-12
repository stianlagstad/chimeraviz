#' Create a reference index using Rsubread.
#'
#' This functions builds an index from a reference sequence using Rsubread.
#'
#' See the plot-fusion-reads vignette for a complete example of how to use the
#' Rsubread functions.
#'
#' @seealso rsubreadAlign
#'
#' @param referenceFasta The path to the reference sequence to build index from.
#'
#' @examples
#' \dontrun{
#' rsubreadIndex(
#'   referenceFasta = "reference.fa")
#' }
#'
#' @importFrom Rsubread buildindex
#'
#' @export
rsubreadIndex <- function(referenceFasta) {

  # Check if the input paths actually exist
  if (file.exists(referenceFasta) == FALSE) {
    stop("Invalid referenceFasta filename. File doesn't exist.")
  }

  Rsubread::buildindex(
    basename = referenceFasta,
    reference = referenceFasta)

  # Did we actually get an output?
  thisFileShouldExist <- paste(referenceFasta, ".00.b.array", sep = "")
  if (file.exists(thisFileShouldExist) == FALSE) {
    warning("It doesn't look like rsubreadIndex() worked.")
  }
}

#' Align reads to fusion sequence using Rsubread.
#'
#' This function aligns paired end reads to the given reference, producing .bam
#' and .bai files. Note that this function also uses Rsamtools to index and sort
#' the .bamfile that Rsubread::align() creates.
#'
#' See the plot-fusion-reads vignette for a complete example of how to use the
#' Rsubread functions.
#'
#' @seealso rsubreadIndex
#'
#' @param referenceName Name of the reference file.
#' @param fastq1 The first fastq file.
#' @param fastq2 The second fastq file.
#' @param outputBamFilename Filename for the resulting .bam file.
#'
#' @examples
#' \dontrun{
#' rsubreadAlign(
#'   referenceName = "reference.fa",
#'   fastq1 = "reads1.fq",
#'   fastq2 = "reads2.fq",
#'   outputBamFilename = "outputbam")
#' }
#'
#' @importFrom Rsamtools sortBam indexBam
#' @importFrom Rsubread align
#'
#' @export
rsubreadAlign <- function(
  referenceName,
  fastq1,
  fastq2,
  outputBamFilename) {

  # Check if the input paths actually exist
  if (file.exists(referenceName) == FALSE) {
    stop("Invalid referenceName filename. File doesn't exist.")
  }
  if (file.exists(fastq1) == FALSE || file.exists(fastq2) == FALSE) {
    stop("Invalid fastq filenames. File(s) doesn't exist.")
  }

  # Align
  Rsubread::align(
    index = referenceName,
    readfile1 = fastq1,
    readfile2 = fastq2,
    output_file = paste(outputBamFilename, ".bam", sep=""))

  # Sort
  Rsamtools::sortBam(
    file = paste(outputBamFilename, ".bam", sep=""),
    destination = outputBamFilename
  )

  # Index
  Rsamtools::indexBam(files = paste(outputBamFilename, ".bam", sep=""))

  # Did we actually get an output?
  if (file.exists(paste(outputBamFilename, ".bam", sep = "")) == FALSE) {
    warning("It doesn't look like rsubreadAlign() worked.")
  }
}
