#' Import results from a JAFFA run into a list of Fusion objects.
#'
#' A function that imports the results from a JAFFA run, typically from a
#' jaffa_results.csv file, into a list of Fusion objects.
#'
#' @param filename Filename for the jaffa_results.csv file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' jaffaData <- system.file(
#'   "extdata",
#'   "jaffa_results.csv",
#'   package = "chimeraviz")
#' fusions <- importJaffa(jaffaData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importJaffa <- function (filename, genomeVersion, limit) {

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
      col_types_jaffa = readr::cols_only(
        "sample" = col_character(),
        "fusion genes" = col_character(),
        "chrom1" = col_character(),
        "base1" = col_integer(),
        "chrom2" = col_character(),
        "base2" = col_integer(),
        "gap (kb)" = col_skip(),
        "spanning pairs" = col_character(),
        "spanning reads" = col_integer(),
        "inframe" = col_character(),
        "aligns" = col_skip(),
        "rearrangement" = col_skip(),
        "contig" = col_skip(),
        "contig break" = col_skip(),
        "classification" = col_character(),
        "known" = col_skip()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_csv(
          file = filename,
          col_types = col_types_jaffa
        )
      } else {
        # Only read up to the limit
        readr::read_csv(
          file = filename,
          col_types = col_types_jaffa,
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
  fusionTool            <- "jaffa"
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import jaffa-specific fields
    fusionToolSpecificData <- list()
    fusionToolSpecificData[["classification"]] <- report[[i, "classification"]]

    # id for this fusion
    id <- as.character(i)

    # Is the downstream fusion partner in-frame?
    if (is.na(report[[i, "inframe"]])) {
      inframe = NA
    }else if(report[[i, "inframe"]]) {
      inframe = TRUE
    } else {
      inframe = FALSE
    }

    # Chromosome names
    chromosomeA <- report[[i, "chrom1"]]
    chromosomeB <- report[[i, "chrom2"]]

    # Breakpoints
    breakpointA <- report[[i, "base1"]]
    breakpointB <- report[[i, "base2"]]

    # Strand
    strandA <- "*"
    strandB <- "*"

    # Number of supporting reads
    splitReadsCount <- as.numeric(report[[i, "spanning reads"]])
    # In the JAFFA example run output the "spanning pairs" field had a couple of
    # "-" values. This causes problems, as we need this to be numeric. If it's
    # indeed a character, just set it to 0:
    spanningReadsCount <- tryCatch(
      spanningReadsCount <- as.numeric(report[[i, "spanning pairs"]]),
      warning = function(w) {
        spanningReadsCount <- 0
        }
    )

    # Not including the fusion junction sequence for JAFFA
    junctionSequenceA <- Biostrings::DNAString()
    junctionSequenceB <- Biostrings::DNAString()

    # Gene names
    geneNames <- unlist(strsplit(report[[i, "fusion genes"]], split=":"))
    nameA <- geneNames[1]
    nameB <- geneNames[2]

    # Ensembl ids
    # ensemblIdA <- geneNames1[2]
    # ensemblIdB <- geneNames2[2]

    # Until we have example data we are sure about, don't set ensemblIds:
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
