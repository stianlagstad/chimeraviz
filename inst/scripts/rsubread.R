#' Create a reference index using Rsubread.
#'
#' This functions builds an index from a reference sequence using Rsubread.
#'
#' See the plot-fusion-reads vignette for a complete example of how to use the
#' Rsubread functions.
#'
#' @seealso rsubread_align
#'
#' @param reference_fasta The path to the reference sequence to build index from.
#'
#' @examples
#' \dontrun{
#' rsubread_index(
#'   reference_fasta = "reference.fa")
#' }
#'
#' @importFrom Rsubread buildindex
#'
#' @export
rsubread_index <- function(reference_fasta) {

  # Check if the input paths actually exist
  if (file.exists(reference_fasta) == FALSE) {
    stop("Invalid referenceFasta filename. File doesn't exist.")
  }

  Rsubread::buildindex(
    basename = reference_fasta,
    reference = reference_fasta)

  # Did we actually get an output?
  this_file_should_exist <- paste(reference_fasta, ".00.b.array", sep = "")
  if (file.exists(this_file_should_exist) == FALSE) {
    warning("It doesn't look like rsubread_index() worked.")
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
#' @seealso rsubread_index
#'
#' @param reference_name Name of the reference file.
#' @param fastq1 The first fastq file.
#' @param fastq2 The second fastq file.
#' @param output_bam_filename Filename for the resulting .bam file.
#'
#' @examples
#' \dontrun{
#' rsubread_align(
#'   reference_name = "reference.fa",
#'   fastq1 = "reads1.fq",
#'   fastq2 = "reads2.fq",
#'   output_bam_filename = "outputbam")
#' }
#'
#' @importFrom Rsamtools sortBam indexBam
#' @importFrom Rsubread align
#'
#' @export
rsubread_align <- function(
  reference_name,
  fastq1,
  fastq2,
  output_bam_filename) {

  # Check if the input paths actually exist
  if (file.exists(reference_name) == FALSE) {
    stop("Invalid referenceName filename. File doesn't exist.")
  }
  if (file.exists(fastq1) == FALSE || file.exists(fastq2) == FALSE) {
    stop("Invalid fastq filenames. File(s) doesn't exist.")
  }

  # Align
  Rsubread::align(
    index = reference_name,
    readfile1 = fastq1,
    readfile2 = fastq2,
    output_file = paste(output_bam_filename, ".bam", sep = ""))

  # Sort
  Rsamtools::sortBam(
    file = paste(output_bam_filename, ".bam", sep = ""),
    destination = output_bam_filename
  )

  # Index
  Rsamtools::indexBam(files = paste(output_bam_filename, ".bam", sep = ""))

  # Did we actually get an output?
  if (file.exists(paste(output_bam_filename, ".bam", sep = "")) == FALSE) {
    warning("It doesn't look like rsubreadAlign() worked.")
  }
}
