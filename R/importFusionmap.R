#' Import results from a FusionMap run into a list of Fusion objects.
#'
#' A function that imports the results from a FusionMap run, typically from
#' a InputFastq.FusionReport.txt file, into a list of Fusion objects.
#'
#' @param filename Filename for the FusionMap PairedEndFusionReport.txt results
#' file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' fusionmapData <- system.file(
#'   "extdata",
#'   "FusionMap_01_TestDataset_InputFastq.FusionReport.txt",
#'   package = "chimeraviz")
#' fusions <- importFusionmap(fusionmapData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importFusionmap <- function (filename, genomeVersion, limit) {

  # Is the genome version valid?
  validGenomes <- c("hg19", "hg38", "mm10")
  if (is.na(match(tolower(genomeVersion), tolower(validGenomes)))) {
    stop("Invalid genome version given")
  }

  # If the limit is set, is the value valid?
  if (missing(limit) == FALSE) {
    if (is.numeric(limit) == FALSE || limit <= 0) {
      stop("limit must be a numeric value bigger than 0")
    }
  }

  # Try to read the fusion report
  report <- withCallingHandlers(
    {
      col_types_fusionmap = readr::cols_only(
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
      message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      warning(cond)
    }
  )

  # Set variables
  id                    <- NA
  inframe               <- NA
  fusionTool            <- "fusionmap"
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import fusionmap-specific fields
    fusionToolSpecificData <- list()
    fusionToolSpecificData[["FrameShift"]] <- report[[i, "FrameShift"]]

    # id for this fusion
    id <- report[[i, "FusionID"]]

    # Is the downstream fusion partner in-frame?
    if(report[[i, "FrameShiftClass"]] == "InFrame") {
      inframe = TRUE
    } else {
      inframe = FALSE
    }

    # Chromosome names
    chromosomeA <- paste("chr", report[[i, "Chromosome1"]], sep = "")
    chromosomeB <- paste("chr", report[[i, "Chromosome2"]], sep = "")

    # Breakpoints
    breakpointA <- report[[i, "Position1"]]
    breakpointB <- report[[i, "Position2"]]

    # Strand
    strands <- unlist(strsplit(report[[i, "Strand"]], split=""))
    strandA <- strands[1]
    strandB <- strands[2]

    # Number of supporting reads. We only get one value from a FusionMap P
    # results file, so store that in splitReadsCount
    splitReadsCount <- report[[i, 2]]
    spanningReadsCount <- 0

    # Get the fusion sequence. Split it into the part from gene1 and gene2
    junctionSequence <- strsplit(report[[i, "FusionJunctionSequence"]], "\\@")
    junctionSequenceA <- Biostrings::DNAString(junctionSequence[[1]][1])
    junctionSequenceB <- Biostrings::DNAString(junctionSequence[[1]][2])

    # Gene names
    nameA <- report[[i, "KnownGene1"]]
    nameB <- report[[i, "KnownGene2"]]

    # Ensembl ids
    ensemblIdA <- NA_character_
    ensemblIdB <- NA_character_

    # PartnerGene objects
    geneA <- new(Class = "PartnerGene",
                 name = nameA,
                 ensemblId = ensemblIdA,
                 chromosome = chromosomeA,
                 breakpoint = breakpointA,
                 strand = strandA,
                 junctionSequence = junctionSequenceA,
                 transcripts = GenomicRanges::GRangesList())
    geneB <- new(Class = "PartnerGene",
                 name = nameB,
                 ensemblId = ensemblIdB,
                 chromosome = chromosomeB,
                 breakpoint = breakpointB,
                 strand = strandB,
                 junctionSequence = junctionSequenceB,
                 transcripts = GenomicRanges::GRangesList())

    fusionList[[i]] <- new(Class = "Fusion",
                           id = id,
                           fusionTool = fusionTool,
                           genomeVersion = genomeVersion,
                           spanningReadsCount = spanningReadsCount,
                           splitReadsCount = splitReadsCount,
                           fusionReadsAlignment = Gviz::AlignmentsTrack(),
                           geneA = geneA,
                           geneB = geneB,
                           inframe = inframe,
                           fusionToolSpecificData = fusionToolSpecificData)
  }

  # Return the list of Fusion objects
  fusionList
}
