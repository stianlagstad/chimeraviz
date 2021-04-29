#' Import results from a SQUID run into a list of Fusion objects.
#'
#' A function that imports the results from a SQUID run into a list of Fusion
#' objects.
#'
#' @param filename Filename for the SQUID results.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' squidfile <- system.file(
#'   "extdata",
#'   "squid_hcc1954_sv.txt",
#'   package="chimeraviz")
#' fusions <- import_squid(squidfile, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.

#' @importFrom data.table fread
#'
#' @export
import_squid <- function(filename, genome_version, limit) {

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
        "chrom1" = "character",
        "start1" = "integer",
        "end1" = "integer",
        "chrom2" = "character",
        "start2" = "integer",
        "end2" = "integer",
        "score" = "integer",
        "strand1" = "character",
        "strand2" = "character",
        "num_concordantfrag_bp1" = "integer",
        "num_concordantfrag_bp2" = "integer"
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
  fusion_tool          <- "SQUID"
  spanning_reads_count <- NA
  split_reads_count    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import defuse-specific fields
    fusion_tool_specific_data <- list()
    fusion_tool_specific_data[["score"]] <- report[[i, "score"]]
    fusion_tool_specific_data[["score"]] <-
      report[[i, "num_concordantfrag_bp1"]]
    fusion_tool_specific_data[["score"]] <-
      report[[i, "num_concordantfrag_bp2"]]

    # Fusion id
    id <- as.character(i)

    # Is the downstream fusion partner in-frame?
    inframe <- NA

    # Strand
    strand_upstream <- report[[i, "strand1"]]
    strand_downstream <- report[[i, "strand2"]]

    # Number of supporting reads
    split_reads_count <- 0
    spanning_reads_count <- 0

    # No junction sequence is given
    junction_sequence_upstream <-
      Biostrings::DNAString()
    junction_sequence_downstream <-
      Biostrings::DNAString()

    # Breakpoints
    if (strand_upstream == "+") {
      breakpoint_upstream <- report[[i, "end1"]]
    } else {
      breakpoint_upstream <- report[[i, "start1"]]
    }
    if (strand_downstream == "+") {
      breakpoint_downstream <- report[[i, "start2"]]
    } else {
      breakpoint_downstream <- report[[i, "end2"]]
    }

    # Chromosome names
    chromosome_upstream <- report[[i, "chrom1"]]
    chromosome_downstream <- report[[i, "chrom2"]]

    # Gene names
    name_upstream <- NA_character_
    name_downstream <- NA_character_

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
