#' Import results from a deFuse run into a list of Fusion objects.
#'
#' A function that imports the results from a deFuse run, typically from a
#' results.filtered.tsv file, into a list of Fusion objects.
#'
#' @param filename Filename for the deFuse results .tsv file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importDefuse <- function (filename, genomeVersion, limit) {

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
      col_types_defuse = readr::cols_only(
        "cluster_id" = readr::col_character(),
        "splitr_sequence" = readr::col_character(),
        "splitr_count" = readr::col_integer(),
        "gene1" = readr::col_character(),
        "gene2" = readr::col_character(),
        "gene_chromosome1" = readr::col_character(),
        "gene_chromosome2" = readr::col_character(),
        "gene_name1" = readr::col_character(),
        "gene_name2" = readr::col_character(),
        "gene_strand1" = readr::col_character(),
        "gene_strand2" = readr::col_character(),
        "genomic_break_pos1" = readr::col_integer(),
        "genomic_break_pos2" = readr::col_integer(),
        "span_count" = readr::col_integer(),
        "orf" = readr::col_character(),
        "probability" = readr::col_number()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_defuse
        )
      } else {
        # Only read up to the limit
        # readr::read_tsv(file = filename, n_max = limit)
        readr::read_tsv(
          file = filename,
          col_types = col_types_defuse,
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
  fusionTool            <- "defuse"
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import defuse-specific fields
    fusionToolSpecificData <- list()
    fusionToolSpecificData[["probability"]] = report[[i, "probability"]]

    # Cluster id
    id <- report[[i, "cluster_id"]]

    # Is the downstream fusion partner in-frame?
    if(report[[i, "orf"]] == "Y") {
      inframe = TRUE
    } else {
      inframe = FALSE
    }

    # Strand
    strandA <- report[[i, "gene_strand1"]]
    strandB <- report[[i, "gene_strand2"]]

    # Number of supporting reads
    splitReadsCount <- report[[i, "splitr_count"]]
    spanningReadsCount <- report[[i, "span_count"]]

    # Get the fusion sequence. Split it into the part from gene1 and gene2
    junctionSequence <- strsplit(report[[i, "splitr_sequence"]], "\\|")
    junctionSequenceA <- Biostrings::DNAString(junctionSequence[[1]][1])
    junctionSequenceB <- Biostrings::DNAString(junctionSequence[[1]][2])

    # Breakpoints
    breakpointA <- report[[i, "genomic_break_pos1"]]
    breakpointB <- report[[i, "genomic_break_pos2"]]

    # Chromosome names
    chromosomeA <- paste("chr", report[[i, "gene_chromosome1"]], sep = "")
    chromosomeB <- paste("chr", report[[i, "gene_chromosome2"]], sep = "")

    # Gene names
    nameA <- report[[i, "gene_name1"]]
    nameB <- report[[i, "gene_name2"]]

    # Ensembl ids
    ensemblIdA <- report[[i, "gene1"]]
    ensemblIdB <- report[[i, "gene2"]]

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
