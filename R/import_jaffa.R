#' Import results from a JAFFA run into a list of Fusion objects.
#'
#' A function that imports the results from a JAFFA run, typically from a
#' jaffa_results.csv file, into a list of Fusion objects.
#'
#' @param filename Filename for the jaffa_results.csv file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' jaffaData <- system.file(
#'   "extdata",
#'   "jaffa_results.csv",
#'   package = "chimeraviz")
#' fusions <- import_jaffa(jaffaData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#' 
#' @importFrom data.table fread
#'
#' @export
import_jaffa <- function (filename, genome_version, limit) {

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
        "sample" = "character",
        "fusion genes" = "character",
        "chrom1" = "character",
        "base1" = "integer",
        "chrom2" = "character",
        "base2" = "integer",
        "spanning pairs" = "character",
        "spanning reads" = "integer",
        "inframe" = "character",
        "classification" = "character"
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
      #message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      # End Exclude Linting
    }
  )

  # Set variables
  id                   <- NA
  inframe              <- NA
  fusion_tool          <- "jaffa"
  spanning_reads_count <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import jaffa-specific fields
    fusion_tool_specific_data <- list()
    fusion_tool_specific_data[["classification"]] <-
      report[[i, "classification"]]

    # id for this fusion
    id <- as.character(i)

    # Is the downstream fusion partner in-frame?
    if (is.na(report[[i, "inframe"]])) {
      inframe <- NA
    } else if (report[[i, "inframe"]]) {
      inframe <- TRUE
    } else {
      inframe <- FALSE
    }

    # Chromosome names
    chromosome_upstream <- report[[i, "chrom1"]]
    chromosome_downstream <- report[[i, "chrom2"]]

    # Breakpoints
    breakpoint_upstream <- report[[i, "base1"]]
    breakpoint_downstream <- report[[i, "base2"]]

    # Strand
    strand_upstream <- "*"
    strand_downstream <- "*"

    # Number of supporting reads
    split_reads_count <- as.numeric(report[[i, "spanning reads"]])
    # In the JAFFA example run output the "spanning pairs" field had a couple of
    # "-" values. This causes problems, as we need this to be numeric. If it's
    # indeed a character, just set it to 0:
    spanning_reads_count <- tryCatch(
      as.numeric(report[[i, "spanning pairs"]]),
      warning = function(w) {
        0
      }
    )

    # Not including the fusion junction sequence for JAFFA
    junction_sequence_upstream <- Biostrings::DNAString()
    junction_sequence_downstream <- Biostrings::DNAString()

    # Gene names
    gene_names <- unlist(strsplit(report[[i, "fusion genes"]], split = ":"))
    name_upstream <- gene_names[1]
    name_downstream <- gene_names[2]

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
