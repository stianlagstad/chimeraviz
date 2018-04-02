#' Import results from a FusionMap run into a list of Fusion objects.
#'
#' A function that imports the results from a FusionMap run, typically from
#' a InputFastq.FusionReport.txt file, into a list of Fusion objects.
#'
#' @param filename Filename for the FusionMap PairedEndFusionReport.txt results
#' file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' fusionmapData <- system.file(
#'   "extdata",
#'   "FusionMap_01_TestDataset_InputFastq.FusionReport.txt",
#'   package = "chimeraviz")
#' fusions <- import_fusionmap(fusionmapData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
import_fusionmap <- function (filename, genome_version, limit) {

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
      col_types_fusionmap <- readr::cols_only(
        "FusionID" = col_character(),
        "DatasetP2_SimulatedReads.UniqueCuttingPositionCount" = col_skip(),
        "DatasetP2_SimulatedReads.SeedCount" = col_skip(),
        "DatasetP2_SimulatedReads.RescuedCount" = col_integer(),
        "Strand" = col_character(),
        "Chromosome1" = col_character(),
        "Position1" = col_integer(),
        "Chromosome2" = col_character(),
        "Position2" = col_integer(),
        "KnownGene1" = col_character(),
        "KnownTranscript1" = col_skip(),
        "KnownExonNumber1" = col_skip(),
        "KnownTranscriptStrand1" = col_skip(),
        "KnownGene2" = col_character(),
        "KnownTranscript2" = col_skip(),
        "KnownExonNumber2" = col_skip(),
        "KnownTranscriptStrand2" = col_skip(),
        "FusionJunctionSequence" = col_character(),
        "FusionGene" = col_skip(),
        "SplicePattern" = col_skip(),
        "SplicePatternClass" = col_skip(),
        "FrameShift" = col_character(),
        "FrameShiftClass" = col_character(),
        "Distance" = col_skip(),
        "OnExonBoundary" = col_skip(),
        "Filter" = col_skip()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_fusionmap
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_fusionmap,
          n_max = limit
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
  fusion_tool          <- "fusionmap"
  spanning_reads_count <- NA
  split_reads_count    <- NA
  junction_sequence    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import fusionmap-specific fields
    fusion_tool_specific_data <- list()
    fusion_tool_specific_data[["FrameShift"]] <- report[[i, "FrameShift"]]

    # id for this fusion
    id <- report[[i, "FusionID"]]

    # Is the downstream fusion partner in-frame?
    if (report[[i, "FrameShiftClass"]] == "InFrame") {
      inframe <- TRUE
    } else {
      inframe <- FALSE
    }

    # Chromosome names
    chromosome_upstream <- paste("chr", report[[i, "Chromosome1"]], sep = "")
    chromosome_downstream <- paste("chr", report[[i, "Chromosome2"]], sep = "")

    # Breakpoints
    breakpoint_upstream <- report[[i, "Position1"]]
    breakpoint_downstream <- report[[i, "Position2"]]

    # Strand
    strands <- unlist(strsplit(report[[i, "Strand"]], split = ""))
    strand_upstream <- strands[1]
    strand_downstream <- strands[2]

    # Number of supporting reads. We only get one value from a FusionMap P
    # results file, so store that in splitReadsCount
    split_reads_count <- report[[i, 2]]
    spanning_reads_count <- 0

    # Get the fusion sequence. Split it into the part from gene1 and gene2
    junction_sequence <- strsplit(report[[i, "FusionJunctionSequence"]], "\\@")
    junction_sequence_upstream <-
      Biostrings::DNAString(junction_sequence[[1]][1])
    junction_sequence_downstream <-
      Biostrings::DNAString(junction_sequence[[1]][2])

    # Gene names
    name_upstream <- report[[i, "KnownGene1"]]
    name_downstream <- report[[i, "KnownGene2"]]

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
