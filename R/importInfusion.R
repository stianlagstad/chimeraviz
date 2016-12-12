#' Import results from an InFusion run into a list of Fusion objects.
#'
#' A function that imports the results from an InFusion run nto a list of Fusion
#' objects.
#'
#' @param filename Filename for the jaffa_results.csv file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' infusionData <- system.file(
#'   "extdata",
#'   "infusion_fusions.txt",
#'   package = "chimeraviz")
#' fusions <- importInfusion(infusionData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importInfusion <- function (filename, genomeVersion, limit) {

  # Is the genome version valid?
  validGenomes <- c("hg19", "hg38")
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
  report <- tryCatch(
    {
      col_types_infusion = readr::cols_only(
        "#id" = col_integer(),
        "ref1" = col_character(),
        "break_pos1" = col_integer(),
        "region1" = col_skip(),
        "ref2" = col_character(),
        "break_pos2" = col_integer(),
        "region2" = col_skip(),
        "num_span" = col_integer(),
        "num_paired" = col_integer(),
        "genes_1" = col_character(),
        "genes_2" = col_character(),
        "fusion_class" = col_skip()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_infusion
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_infusion,
          n_max = limit
        )
      }
    },
    error = function(cond) {
      message("Reading filename caused an error:")
      stop(cond)
    },
    warning = function(cond) {
      message("Reading filename caused a warning:")
      warning(cond)
    }
  )

  # Set variables
  id                    <- NA
  inframe               <- NA
  fusionTool            <- "infusion"
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import infusion-specific fields
    fusionToolSpecificData <- list()

    # id for this fusion
    id <- as.character(report[[i, "#id"]])

    # Chromosome names
    chromosomeA <- paste("chr", report[[i, "ref1"]], sep = "")
    chromosomeB <- paste("chr", report[[i, "ref2"]], sep = "")

    # Breakpoints
    breakpointA <- report[[i, "break_pos1"]]
    breakpointB <- report[[i, "break_pos2"]]

    # Strand
    strandA <- "*"
    strandB <- "*"

    # Number of supporting reads
    splitReadsCount <- as.numeric(report[[i, "num_span"]])
    spanningReadsCount <- as.numeric(report[[i, "num_paired"]])
    # Note: the terminology used for spanning and split reads varies. According
    # to the source code of InFusion (the file fusion.py), the "num_paired"
    # field is what the author calls "bridge reads", which is what we are
    # calling spanning reads.

    # InFusiondoesn't provide the fusion sequence
    junctionSequenceA <- Biostrings::DNAString()
    junctionSequenceB <- Biostrings::DNAString()

    # Gene names
    nameA <- report[[i, "genes_1"]]
    nameB <- report[[i, "genes_2"]]

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
