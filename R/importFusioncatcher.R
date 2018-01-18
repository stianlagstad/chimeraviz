#' Import results from a Fusioncatcher run into a list of Fusion objects.
#'
#' A function that imports the results from a Fusioncatcher run, typically from
#' a final-list-candidate-fusion-genes.txt file, into a list of Fusion objects.
#'
#' @param filename Filename for the Fusioncatcher
#' final-list-candidate-fusion-genes.txt results file.
#' @param genomeVersion Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' fusioncatcher833ke <- system.file(
#'   "extdata",
#'   "fusioncatcher_833ke_final-list-candidate-fusion-genes.txt",
#'   package = "chimeraviz")
#' fusions <- importFusioncatcher(fusioncatcher833ke, "hg38", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
importFusioncatcher <- function (filename, genomeVersion, limit) {

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
      col_types_fusioncatcher = readr::cols_only(
        "Gene_1_symbol(5end_fusion_partner)" = readr::col_character(),
        "Gene_2_symbol(3end_fusion_partner)" = readr::col_character(),
        "Spanning_pairs" = readr::col_integer(),
        "Spanning_unique_reads" = readr::col_integer(),
        "Fusion_point_for_gene_1(5end_fusion_partner)" = readr::col_character(),
        "Fusion_point_for_gene_2(3end_fusion_partner)" = readr::col_character(),
        "Gene_1_id(5end_fusion_partner)" = readr::col_character(),
        "Gene_2_id(3end_fusion_partner)" = readr::col_character(),
        "Fusion_sequence" = readr::col_character(),
        "Predicted_effect" = readr::col_character()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_fusioncatcher
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_fusioncatcher,
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
  fusionTool            <- "fusioncatcher"
  spanningReadsCount    <- NA
  splitReadsCount       <- NA
  junctionSequence      <- NA
  fusionReadsAlignment  <- NA

  # List to hold all Fusion objects
  fusionList <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import fusioncatcher-specific fields
    fusionToolSpecificData <- list()

    # Set id for this fusion. Fusioncatcher does not have any unique identifier
    # for each row of the results file, so just use our counter.
    id <- as.character(i)

    # Is the downstream fusion partner in-frame?
    if(report[[i, "Predicted_effect"]] == "in-frame") {
      inframe = TRUE
    } else {
      inframe = FALSE
    }

    geneAposition <- unlist(strsplit(report[[i, "Fusion_point_for_gene_1(5end_fusion_partner)"]], split=":"))
    geneBposition <- unlist(strsplit(report[[i, "Fusion_point_for_gene_2(3end_fusion_partner)"]], split=":"))

    # Chromosome names
    chromosomeA <- paste("chr", geneAposition[1], sep = "")
    chromosomeB <- paste("chr", geneBposition[1], sep = "")

    # Breakpoints
    breakpointA <- as.numeric(geneAposition[2])
    breakpointB <- as.numeric(geneBposition[2])

    # Strand
    strandA <- geneAposition[3]
    strandB <- geneBposition[3]

    # Number of supporting reads
    splitReadsCount <- report[[i, "Spanning_pairs"]]
    spanningReadsCount <- report[[i, "Spanning_unique_reads"]]

    # Get the fusion sequence. Split it into the part from gene1 and gene2
    junctionSequence <- strsplit(report[[i, "Fusion_sequence"]], split="\\*")
    junctionSequenceA <- Biostrings::DNAString(junctionSequence[[1]][1])
    junctionSequenceB <- Biostrings::DNAString(junctionSequence[[1]][2])

    # Gene names
    nameA <- report[[i, "Gene_1_symbol(5end_fusion_partner)"]]
    nameB <- report[[i, "Gene_2_symbol(3end_fusion_partner)"]]

    # Ensembl ids
    ensemblIdA <- report[[i, "Gene_1_id(5end_fusion_partner)"]]
    ensemblIdB <- report[[i, "Gene_2_id(3end_fusion_partner)"]]

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
