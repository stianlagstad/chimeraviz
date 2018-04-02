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
#' @seealso bowtie_align
#'
#' @param bowtie_build_location The path to bowtie-build on your system.
#' @param reference_fasta The path to the reference sequence to build index from.
#'
#' @examples
#' \dontrun{
#' bowtie_index(
#'   bowtie_build_location = "path/to/bowtie-build",
#'   reference_fasta = "reference.fa")
#' }
#'
#' @export
bowtie_index <- function(bowtie_build_location, reference_fasta) {

  # Check if the input paths actually exist
  if (file.exists(bowtie_build_location) == FALSE) {
    stop("Invalid bowtie_build_location File doesn't exist.")
  }
  if (file.exists(reference_fasta) == FALSE) {
    stop("Invalid reference_fasta filename. File doesn't exist.")
  }

  create_index_command <- paste(
    shQuote(bowtie_build_location),
    shQuote(reference_fasta),
    shQuote(reference_fasta),
    sep = " "
  )

  # This should create a command like this:
  # "bowtie2-build referenceFasta reference_fasta" # Exclude Linting

  # From the bowtie manual: "The command should print many lines of output then
  # quit. When the command completes, the current directory will contain four
  # new files that all start with `reference_fasta` and end with `.1.bt2`,
  # `.2.bt2`, `.3.bt2`, `.4.bt2`, `.rev.1.bt2`, and `.rev.2.bt2`. These files
  # constitute the index."

  # Run the command and wait for it to complete
  system(create_index_command, wait = TRUE)

  # Did we actually get an output?
  this_file_should_exist <- paste(reference_fasta, ".1.bt2", sep = "")
  if (file.exists(this_file_should_exist) == FALSE) {
    warning(paste("It doesn't look like bowtie_index() worked. Command tried:",
                  create_index_command))
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
#' @seealso bowtie_index
#'
#' @param bowtie_location The path to bowtie on your system.
#' @param reference_name Name of the reference file.
#' @param fastq1 The first fastq file.
#' @param fastq2 The second fastq file.
#' @param output_bam_filename Filename for the resulting .bam file.
#' @param p Number of alignment threads to launch (default: 1). Given to -p (or
#' --threads).
#'
#' @examples
#' \dontrun{
#' bowtie_align(
#'   bowtie_location = "path/to/bowtie",
#'   reference_name = "reference.fa",
#'   fastq1 = "reads1.fq",
#'   fastq2 = "reads2.fq",
#'   output_bam_filename = "outputbam",
#'   p = 2)
#' }
#'
#' @export
bowtie_align <- function(
  bowtie_location,
  reference_name,
  fastq1,
  fastq2,
  output_bam_filename,
  p = 1
) {

  # Check if the input paths actually exist
  if (file.exists(bowtie_location) == FALSE) {
    stop("Invalid bowtie_location File doesn't exist.")
  }
  if (file.exists(reference_name) == FALSE) {
    stop("Invalid reference_name filename. File doesn't exist.")
  }
  if (file.exists(fastq1) == FALSE || file.exists(fastq2) == FALSE) {
    stop("Invalid fastq filenames. File(s) doesn't exist.")
  }
  if (!any(class(p) == c("integer", "numeric"))) {
    stop("Invalid p. Must be integer.")
  }

  align_command <- paste(
    bowtie_location,
    "-x", shQuote(reference_name),
    "-1", shQuote(fastq1),
    "-2", shQuote(fastq2),
    "--very-sensitive", "-p", p,
    "-S",
    paste(output_bam_filename, ".sam", sep = ""),
    sep = " "
  )

  # This should create a command like this:
  # "bowtie2 -x ref.fa -1 fastq1 -2 fastq2 --very-sensitive -p 1

  # Run the command and wait for it to complete
  system(align_command, wait = TRUE)

  # Sort and index
  Rsamtools::asBam(
    file = paste(output_bam_filename, ".sam", sep = ""),
    destination = output_bam_filename)

  # Did we actually get an output?
  if (file.exists(paste(output_bam_filename, ".bam", sep = "")) == FALSE) {
    warning(paste("It doesn't look like bowtie_align() worked. Command tried:",
                  align_command))
  }
}
