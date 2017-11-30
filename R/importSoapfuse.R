#' Import results from a SOAPfuse run into a list of Fusion objects.
#'
#' A function that imports the results from a SOAPfuse run, typically from
#' a final.Fusion.specific.for.genes file, into a list of Fusion objects.
#'
#' @param filename Filename for the SOAPfuse
#' final-list-candidate-fusion-genes.txt results file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' soapfuse833ke <- system.file(
#'   "extdata",
#'   "soapfuse_833ke_final.Fusion.specific.for.genes",
#'   package = "chimeraviz")
#' fusions <- importSoapfuse(soapfuse833ke, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importSoapfuse <- function (filename, genomeVersion, limit) {

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
      col_types_soapfuse = readr::cols_only(
        "up_gene" = readr::col_character(),
        "up_chr" = readr::col_character(),
        "up_strand" = readr::col_character(),
        "up_Genome_pos" = readr::col_integer(),
        "dw_gene" = readr::col_character(),
        "dw_chr" = readr::col_character(),
        "dw_strand" = readr::col_character(),
        "dw_Genome_pos" = readr::col_integer(),
        "Span_reads_num" = readr::col_integer(),
        "Junc_reads_num" = readr::col_integer(),
        "down_fusion_part_frame-shift_or_not" = readr::col_character()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_soapfuse
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_soapfuse,
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
  fusionTool            <- "soapfuse"
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import soapfuse-specific fields
    fusionToolSpecificData <- list()

    # Set id for this fusion. Soapfuse does not have any unique identifier
    # for each row of the results file, so just use our counter.
    id <- as.character(i)

    # Is the downstream fusion partner in-frame?
    if (is.na(report[[i, "down_fusion_part_frame-shift_or_not"]])) {
      inframe = NA
    }else if(report[[i, "down_fusion_part_frame-shift_or_not"]] == "inframe-shift") {
      inframe = TRUE
    } else {
      inframe = FALSE
    }

    # Chromosome names
    chromosomeA <- report[[i, "up_chr"]]
    chromosomeB <- report[[i, "dw_chr"]]

    # Breakpoints
    breakpointA <- report[[i, "up_Genome_pos"]]
    breakpointB <- report[[i, "dw_Genome_pos"]]

    # Strand
    strandA <- report[[i, "up_strand"]]
    strandB <- report[[i, "dw_strand"]]

    # Number of supporting reads
    splitReadsCount <- report[[i, "Span_reads_num"]]
    spanningReadsCount <- report[[i, "Junc_reads_num"]]

    # SOAPfuse doesn't provide the fusion sequence
    junctionSequenceA <- Biostrings::DNAString()
    junctionSequenceB <- Biostrings::DNAString()

    # Gene names
    nameA <- report[[i, "up_gene"]]
    nameB <- report[[i, "dw_gene"]]

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
