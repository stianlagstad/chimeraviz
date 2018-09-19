#' Import results from an InFusion run into a list of Fusion objects.
#'
#' A function that imports the results from an InFusion run nto a list of Fusion
#' objects.
#'
#' @param filename Filename for the jaffa_results.csv file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' infusionData <- system.file(
#'   "extdata",
#'   "infusion_fusions.txt",
#'   package = "chimeraviz")
#' fusions <- import_infusion(infusionData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @importFrom data.table fread
#'
#' @export
import_infusion <- function (filename, genome_version, limit) {

  # Is the genome version valid?
  valid_genomes <- c("hg19", "hg38", "mm10")
  if (is.na(match(tolower(genome_version), tolower(valid_genomes)))) {
    stop("Invalid genome version given")
  }

  # If the limit is set, is the value valid?
  if (missing(limit) == FALSE) {
    if (is.numeric(limit) == FALSE || limit <= 0) {
      stop("limit must be a numeric value bigger than 0")
    }
  }

  # Try to read the fusion report
  report <- withCallingHandlers({
      col_types <- c(
        "#id" = "integer",
        "ref1" = "character",
        "break_pos1" = "integer",
        "ref2" = "character",
        "break_pos2" = "integer",
        "num_span" = "integer",
        "num_paired" = "integer",
        "genes_1" = "character",
        "genes_2" = "character"
      )
      if (missing(limit)) {
        # Read all lines
        data.table::fread(
          input = filename,
          colClasses = col_types,
          showProgress = FALSE
        )
      } else {
        # Only read up to the limit
        data.table::fread(
          input = filename,
          colClasses = col_types,
          showProgress = FALSE,
          nrows = limit
        )
      }
    },
    error = function(cond) {
      message(paste0("Reading ", filename, " caused an error: ", cond[[1]]))
      stop(cond)
    },
    warning = function(cond) {
      # Begin Exclude Linting
      message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      # End Exclude Linting
    }
  )

  # Set variables
  id                   <- NA
  inframe              <- NA
  fusion_tool          <- "infusion"
  spanning_reads_count <- NA
  split_reads_count    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import infusion-specific fields
    fusion_tool_specific_data <- list()

    # id for this fusion
    id <- as.character(report[[i, "#id"]])

    # Chromosome names
    chromosome_upstream <- paste("chr", report[[i, "ref1"]], sep = "")
    chromosome_downstream <- paste("chr", report[[i, "ref2"]], sep = "")

    # Breakpoints
    breakpoint_upstream <- report[[i, "break_pos1"]]
    breakpoint_downstream <- report[[i, "break_pos2"]]

    # Strand
    strand_upstream <- "*"
    strand_downstream <- "*"

    # Number of supporting reads
    split_reads_count <- as.numeric(report[[i, "num_span"]])
    spanning_reads_count <- as.numeric(report[[i, "num_paired"]])
    # Note: the terminology used for spanning and split reads varies. According
    # to the source code of InFusion (the file fusion.py), the "num_paired"
    # field is what the author calls "bridge reads", which is what we are
    # calling spanning reads.

    # InFusiondoesn't provide the fusion sequence
    junction_sequence_upstream <- Biostrings::DNAString()
    junction_sequence_downstream <- Biostrings::DNAString()

    # Gene names
    name_upstream <- report[[i, "genes_1"]]
    name_downstream <- report[[i, "genes_2"]]

    # Ensembl ids
    ensembl_id_upstream <- NA_character_
    ensembl_id_downstream <- NA_character_

    # PartnerGene objects
    gene_upstream <- new(
      Class = "PartnerGene",
      name = name_upstream,
      ensembl_id = ensembl_id_upstream,
      chromosome = chromosome_upstream,
      breakpoint = breakpoint_upstream,
      strand = strand_upstream,
      junction_sequence = junction_sequence_upstream,
      transcripts = GenomicRanges::GRangesList()
    )
    gene_downstream <- new(
      Class = "PartnerGene",
      name = name_downstream,
      ensembl_id = ensembl_id_downstream,
      chromosome = chromosome_downstream,
      breakpoint = breakpoint_downstream,
      strand = strand_downstream,
      junction_sequence = junction_sequence_downstream,
      transcripts = GenomicRanges::GRangesList()
    )

    fusion_list[[i]] <- new(
      Class = "Fusion",
      id = id,
      fusion_tool = fusion_tool,
      genome_version = genome_version,
      spanning_reads_count = spanning_reads_count,
      split_reads_count = split_reads_count,
      fusion_reads_alignment = Gviz::AlignmentsTrack(),
      gene_upstream = gene_upstream,
      gene_downstream = gene_downstream,
      inframe = inframe,
      fusion_tool_specific_data = fusion_tool_specific_data
    )
  }

  # Return the list of Fusion objects
  fusion_list
}
