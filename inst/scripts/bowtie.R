#' Create a reference index using bowtie-build
#'
#' This functions builds an index from a reference sequence using bowtie.
#'
#' Both bowtie and bowtie2 can be used.
#'
#' From the bowtie manual: "The command should print many lines of output then
#' quit. When the command completes, the current directory will contain four new
#' files that all start with `referenceFasta` and end with `.1.bt2`, `.2.bt2`,
#' `.3.bt2`, `.4.bt2`, `.rev.1.bt2`, and `.rev.2.bt2`. These files constitute
#' the index."
#'
#' See the plot-fusion-reads vignette for a complete example of how to use the
#' Bowtie functions.
#'
#' @seealso bowtieAlign
#'
#' @param bowtieBuildLocation The path to bowtie-build on your system.
#' @param referenceFasta The path to the reference sequence to build index from.
#'
#' @examples
#' \dontrun{
#' bowtieIndex(
#'   bowtieBuildLocation = "path/to/bowtie-build",
#'   referenceFasta = "reference.fa")
#' }
#'
#' @export
bowtieIndex <- function(bowtieBuildLocation, referenceFasta) {

  # Check if the input paths actually exist
  if (file.exists(bowtieBuildLocation) == FALSE) {
    stop("Invalid bowtieBuildLocation. File doesn't exist.")
  }
  if (file.exists(referenceFasta) == FALSE) {
    stop("Invalid referenceFasta filename. File doesn't exist.")
  }

  createIndexCommand <- paste(shQuote(bowtieBuildLocation),
                              shQuote(referenceFasta),
                              shQuote(referenceFasta),
                              sep = " ")

  # This should create a command like this:
  # "bowtie2-build referenceFasta referenceFasta"

  # From the bowtie manual: "The command should print many lines of output then
  # quit. When the command completes, the current directory will contain four
  # new files that all start with `referenceFasta` and end with `.1.bt2`,
  # `.2.bt2`, `.3.bt2`, `.4.bt2`, `.rev.1.bt2`, and `.rev.2.bt2`. These files
  # constitute the index."

  # Run the command and wait for it to complete
  system(createIndexCommand, wait = TRUE)

  # Did we actually get an output?
  thisFileShouldExist <- paste(referenceFasta, ".1.bt2", sep = "")
  if (file.exists(thisFileShouldExist) == FALSE) {
    warning(paste("It doesn't look like bowtieIndex() worked. Command tried:",
                  createIndexCommand))
  }
}

#' Align reads to fusion sequence using bowtie
#'
#' This function aligns paired end reads to the given reference, producing .bam
#' and .bai files. Note that this function also uses Rsamtools to index and sort
#' the .sam file that bowtie creates.
#'
#' Both bowtie and bowtie2 can be used.
#'
#' See the plot-fusion-reads vignette for a complete example of how to use the
#' Bowtie functions.
#'
#' @seealso bowtieIndex
#'
#' @param bowtieLocation The path to bowtie on your system.
#' @param referenceName Name of the reference file.
#' @param fastq1 The first fastq file.
#' @param fastq2 The second fastq file.
#' @param outputBamFilename Filename for the resulting .bam file.
#' @param p Number of alignment threads to launch (default: 1). Given to -p (or
#' --threads).
#'
#' @examples
#' \dontrun{
#' bowtieAlign(
#'   bowtieLocation = "path/to/bowtie",
#'   referenceName = "reference.fa",
#'   fastq1 = "reads1.fq",
#'   fastq2 = "reads2.fq",
#'   outputBamFilename = "outputbam",
#'   p = 2)
#' }
#'
#' @export
bowtieAlign <- function(bowtieLocation,
                        referenceName,
                        fastq1,
                        fastq2,
                        outputBamFilename,
                        p = 1) {

  # Check if the input paths actually exist
  if (file.exists(bowtieLocation) == FALSE) {
    stop("Invalid bowtieLocation File doesn't exist.")
  }
  if (file.exists(referenceName) == FALSE) {
    stop("Invalid referenceName filename. File doesn't exist.")
  }
  if (file.exists(fastq1) == FALSE || file.exists(fastq2) == FALSE) {
    stop("Invalid fastq filenames. File(s) doesn't exist.")
  }
  if (!any(class(p) == c("integer", "numeric"))) {
    stop("Invalid p. Must be integer.")
  }

  alignCommand <- paste(bowtieLocation,
                        "-x", shQuote(referenceName),
                        "-1", shQuote(fastq1),
                        "-2", shQuote(fastq2),
                        "--very-sensitive", "-p", p,
                        "-S",
                        paste(outputBamFilename, ".sam", sep = ""),
                        sep = " ")

  # This should create a command like this:
  # "bowtie2 -x ref.fa -1 fastq1 -2 fastq2 --very-sensitive -p 1

  # Run the command and wait for it to complete
  system(alignCommand, wait = TRUE)

  # Sort and index
  Rsamtools::asBam(
    file = paste(outputBamFilename, ".sam", sep = ""),
    destination = outputBamFilename)

  # Did we actually get an output?
  if (file.exists(paste(outputBamFilename, ".bam", sep = "")) == FALSE) {
    warning(paste("It doesn't look like bowtieAlign() worked. Command tried:",
                  alignCommand))
  }
}
