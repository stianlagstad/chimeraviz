#' Import results from a PRADA run into a list of Fusion objects.
#'
#' A function that imports the results from a PRADA run into a list of Fusion
#' objects.
#'
#' @param filename Filename for the PRADA results file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' pradaData <- system.file(
#'   "extdata",
#'   "PRADA.acc.fusion.fq.TAF.tsv",
#'   package = "chimeraviz")
#' fusions <- importPrada(pradaData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importPrada <- function (filename, genomeVersion, limit) {

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
  report <- withCallingHandlers(
    {
      col_types_prada = readr::cols_only(
        "sample" = col_character(),
        "Gene_A" = col_skip(),
        "Gene_B" = col_skip(),
        "A_chr" = col_skip(),
        "B_chr" = col_skip(),
        "A_strand" = col_integer(),
        "B_strand" = col_integer(),
        "Discordant_n" = col_integer(),
        "JSR_n" = col_integer(),
        "perfectJSR_n" = col_skip(),
        "Junc_n" = col_skip(),
        "Position_Consist" = col_skip(),
        "Junction" = col_character(),
        "Identity" = col_skip(),
        "Align_Len" = col_skip(),
        "Evalue" = col_skip(),
        "BitScore" = col_skip(),
        "TAF_info" = col_skip(),
        "max_TAF" = col_skip()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_prada
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_prada,
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
  fusionTool            <- "prada"
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import prada-specific fields
    fusionToolSpecificData <- list()

    # id for this fusion
    id <- as.character(i)

    # The format of report[[i, "Junction"]] should be:
    # "RBM14:11:66386032_BBS1:11:66283011,6". Split into two first:
    junctionDataGeneA <- strsplit(report[[i, "Junction"]], split="_")[[1]][1]
    junctionDataGeneB <- strsplit(report[[i, "Junction"]], split="_")[[1]][2]
    # The format of junctionDataGeneB should now be: "BBS1:11:66283011,6". So
    # remove everything after comma
    junctionDataGeneB <- substr(
      junctionDataGeneB,
      1,
      unlist(gregexpr(pattern =",", junctionDataGeneB))-1)
    # Now split both on :
    junctionDataGeneA <- unlist(strsplit(junctionDataGeneA, split=":"))
    junctionDataGeneB <- unlist(strsplit(junctionDataGeneB, split=":"))

    # Chromosome names
    chromosomeA <- paste("chr", junctionDataGeneA[[2]], sep = "")
    chromosomeB <- paste("chr", junctionDataGeneB[[2]], sep = "")

    # Breakpoints
    breakpointA <- as.numeric(junctionDataGeneA[[3]])
    breakpointB <- as.numeric(junctionDataGeneB[[3]])

    # Strand
    strandA <- if (report[[i, "A_strand"]] == 1) "+" else "-"
    strandB <- if (report[[i, "B_strand"]] == 1) "+" else "-"

    # Number of supporting reads
    splitReadsCount <- as.numeric(report[[i, "JSR_n"]])
    spanningReadsCount <- as.numeric(report[[i, "Discordant_n"]])

    # PRADA doesn't provide the fusion sequence
    junctionSequenceA <- Biostrings::DNAString()
    junctionSequenceB <- Biostrings::DNAString()

    # Gene names
    nameA <- junctionDataGeneA[[1]]
    nameB <- junctionDataGeneB[[1]]

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
