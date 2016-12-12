#' Import results from a EricScript run into a list of Fusion objects.
#'
#' A function that imports the results from a EricScript run into a list of
#' Fusion objects.
#'
#' @param filename Filename for the EricScript results file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' ericscriptData <- system.file(
#'   "extdata",
#'   "ericscript_SRR1657556.results.total.tsv",
#'   package = "chimeraviz")
#' fusions <- importEricscript(ericscriptData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importEricscript <- function (filename, genomeVersion, limit) {

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
      col_types_ericscript = readr::cols_only(
        "GeneName1" = col_character(),
        "GeneName2" = col_character(),
        "chr1" = col_character(),
        "Breakpoint1" = col_character(),
        "strand1" = col_character(),
        "chr2" = col_character(),
        "Breakpoint2" = col_character(),
        "strand2" = col_character(),
        "EnsemblGene1" = col_character(),
        "EnsemblGene2" = col_character(),
        "crossingreads" = col_integer(),
        "spanningreads" = col_integer(),
        "mean.insertsize" = col_skip(),
        "homology" = col_skip(),
        "fusiontype" = col_skip(),
        "Blacklist" = col_skip(),
        "InfoGene1" = col_skip(),
        "InfoGene2" = col_skip(),
        "JunctionSequence" = col_character(),
        "GeneExpr1" = col_skip(),
        "GeneExpr2" = col_skip(),
        "GeneExpr_Fused" = col_skip(),
        "ES" = col_skip(),
        "GJS" = col_skip(),
        "US" = col_skip(),
        "EricScore" = col_number()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_ericscript
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_ericscript,
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
  fusionTool            <- "ericscript"
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import ericscript-specific fields
    fusionToolSpecificData <- list()
    fusionToolSpecificData[["EricScore"]] <- report[[i, "EricScore"]]

    # Set id for this fusion. Soapfuse does not have any unique identifier
    # for each row of the results file, so just use our counter.
    id <- as.character(i)

    # Chromosome names
    chromosomeA <- paste("chr", report[[i, "chr1"]], sep = "")
    chromosomeB <- paste("chr", report[[i, "chr2"]], sep = "")

    # Breakpoints
    # EricScript doesn't always report fusion breakpoints. If missing, then set
    # breakpoint to -1
    breakpointA <- tryCatch(
      breakpointA <- as.numeric(report[[i, "Breakpoint1"]]),
      warning = function(w) {
        breakpointA <- -1
      }
    )
    breakpointB <- tryCatch(
      breakpointB <- as.numeric(report[[i, "Breakpoint2"]]),
      warning = function(w) {
        breakpointB <- -1
      }
    )

    # Strand
    strandA <- report[[i, "strand1"]]
    strandB <- report[[i, "strand2"]]

    # Number of supporting reads
    splitReadsCount <- report[[i, "crossingreads"]]
    spanningReadsCount <- report[[i, "spanningreads"]]

    # Get the fusion sequence
    # EricScript outputs this as sequencebeforejunctionSEQUENCEAFTERJUNCTION, so
    # first locate the first upper cased letter:
    splitAt <- regexpr(pattern ='([[:upper:]])', report[[i, "JunctionSequence"]])[1]
    # Then we can extract each piece:
    junctionSequenceA <- Biostrings::DNAString(
      substr(
        report[[i, "JunctionSequence"]],
        1,
        splitAt-1)
    )
    junctionSequenceB <- Biostrings::DNAString(
      substr(
        report[[i, "JunctionSequence"]],
        splitAt,
        nchar(report[[i, "JunctionSequence"]]))
    )

    # Gene names
    nameA <- report[[i, "GeneName1"]]
    nameB <- report[[i, "GeneName2"]]

    # Ensembl ids
    ensemblIdA <- report[[i, "EnsemblGene1"]]
    ensemblIdB <- report[[i, "EnsemblGene2"]]

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
