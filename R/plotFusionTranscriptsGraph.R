#' Graph plot of possible fusion transcripts.
#'
#' This function takes a fusion object and a TranscriptDb object and plots a
#' graph showing the possible fusion transcripts.
#'
#' @param fusion The Fusion object to plot.
#' @param edb The edb object that will be used to fetch data.
#' @param whichTranscripts This character vector decides which transcripts are
#' to be plotted. Can be "exonBoundary", "withinExon", "withinIntron",
#' "intergenic", or a character vector with specific transcript ids. Default
#' value is "exonBoundary".
#' @param rankdir Choose whether the graph should be plotted from left to right
#' ("LR"), or from top to bottom ("TB"). This parameter is given to
#' Rgraphviz::plot().
#'
#' @return Creates a fusion transcripts graph plot.
#'
#' @examples
#' # Load data and example fusion event
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 1)
#' fusion <- getFusionById(fusions, 5267)
#' # Load edb
#' edbSqliteFile <- system.file(
#'   "extdata",
#'   "Homo_sapiens.GRCh37.74.sqlite",
#'   package="chimeraviz")
#' edb <- ensembldb::EnsDb(edbSqliteFile)
#' # Temporary file to store the plot
#' pngFilename <- tempfile(
#'   pattern = "fusionPlot",
#'   fileext = ".png",
#'   tmpdir = tempdir())
#' # Open device
#' png(pngFilename, width = 500, height = 500)
#' # Plot!
#' plotFusionTranscriptsGraph(
#'   fusion = fusion,
#'   edb = edb)
#' # Close device
#' dev.off()
#'
#' @export
plotFusionTranscriptsGraph <- function(
  fusion,
  edb = NULL,
  whichTranscripts = "exonBoundary",
  rankdir = "TB") {

  .validatePlotFusionTranscriptsGraphParams(fusion, edb, whichTranscripts, rankdir)
  fusion <- .getTranscriptsIfNotThere(fusion, edb)

  # Select which transcripts to use
  transcriptsA <- selectTranscript(fusion@geneA, whichTranscripts)
  transcriptsB <- selectTranscript(fusion@geneB, whichTranscripts)

  # Unlist
  transcriptsA <- unlist(transcriptsA)
  transcriptsB <- unlist(transcriptsB)

  # calculate fusion transcripts

  # Select the exons that are part of a possible fusion transcript. I.e. these:
  #
  # 1: geneA is on the + strand
  #    Exons with start/end positions before breakpoint
  # 2: geneA is on the - strand
  #    Exons with start/end positions after breakpoint
  # 3: geneB is on the + strand
  #    Exons with start/end positions after breakpoint
  # 4: geneB is on the - strand
  #    Exons with start/end positions before breakpoint

  if (fusion@geneA@strand == "+") {
    fusionExonsA <- transcriptsA[start(ranges(transcriptsA)) < fusion@geneA@breakpoint & end(ranges(transcriptsA)) <= fusion@geneA@breakpoint]
  } else {
    fusionExonsA <- transcriptsA[start(ranges(transcriptsA)) >= fusion@geneA@breakpoint & end(ranges(transcriptsA)) > fusion@geneA@breakpoint]
  }
  if (fusion@geneB@strand == "+") {
    fusionExonsB <- transcriptsB[start(ranges(transcriptsB)) >= fusion@geneB@breakpoint & end(ranges(transcriptsB)) > fusion@geneB@breakpoint]
  } else {
    fusionExonsB <- transcriptsB[start(ranges(transcriptsB)) < fusion@geneB@breakpoint & end(ranges(transcriptsB)) <= fusion@geneB@breakpoint]
  }

  # If either fusionExonsA or fusionExonsB is empty, then the breakpoint does
  # not occur inside any of the transcripts.
  if (length(fusionExonsA) == 0) {
    stop(paste(
      "None of the transcripts given for gene A has the fusion breakpoint ",
      "within them. This plot cannot be created with the given transcripts."))
  }
  if (length(fusionExonsB) == 0) {
    stop(paste(
      "None of the transcripts given for gene B has the fusion breakpoint ",
      "within them. This plot cannot be created with the given transcripts."))
  }

  # For each transcript in fusionExonsA, match the exons in those transcripts ..

  transcriptNamesA <- unique(unlist(mcols(fusionExonsA)$transcript))
  transcriptNamesB <- unique(unlist(mcols(fusionExonsB)$transcript))
  fusionTranscripts <- GRangesList()
  for (i in 1:length(transcriptNamesA)) {
    trNameA <- transcriptNamesA[[i]]
    for (j in 1:length(transcriptNamesB)) {
      trNameB <- transcriptNamesB[[j]]

      exonsA <- fusionExonsA[mcols(fusionExonsA)$transcript == trNameA]
      exonsB <- fusionExonsB[mcols(fusionExonsB)$transcript == trNameB]

      # Reverse order of exons if on minus strand
      if (fusion@geneA@strand == "-") {
        exonsA <- rev(exonsA)
      }
      if (fusion@geneB@strand == "-") {
        exonsB <- rev(exonsB)
      }

      new <- append(
        exonsA,
        exonsB
      )
      mcols(new)$transcript <- paste("fusion", i, j, sep = "")
      fusionTranscripts <- append(
        fusionTranscripts,
        GRangesList(new)
      )
    }
  }

  fusionTranscripts <- unlist(fusionTranscripts)

  exonnames <- unique(mcols(fusionTranscripts)$exon_id)
  exonnamesA <- unique(mcols(fusionTranscripts[mcols(fusionTranscripts)$gene_id == fusion@geneA@ensemblId])$exon_id)
  exonnamesB <- unique(mcols(fusionTranscripts[mcols(fusionTranscripts)$gene_id == fusion@geneB@ensemblId])$exon_id)

  # These lists will hold edge information
  from <- c()
  to <- c()

  transcriptNames <- sort(unique(mcols(fusionTranscripts)$transcript))
  for (i in 1:length(transcriptNames)) {
    exonsInTranscript <- fusionTranscripts[mcols(fusionTranscripts)$transcript == transcriptNames[i]]
    for (j in 2:length(exonsInTranscript)) {
      # if (j == 2) {
      #   # Add start edge
      #   from <- c(from, "start")
      #   to <- c(to, mcols(exonsInTranscript[j-1])$exon_id)
      # }
      # Add normal edge
      from <- c(from, mcols(exonsInTranscript[j-1])$exon_id)
      to <- c(to, mcols(exonsInTranscript[j])$exon_id)
    }
    # # Add stop edge
    # from <- c(from, mcols(exonsInTranscript[j])$exon_id)
    # to <- c(to, "stop")
  }

  # Use the edge information we've gathered and create a graph
  ft <- cbind(from, to)
  # Remove edges that go to the same node. Do this because the way
  # getTranscriptsEnsembldb is implemented it will split exons in which the cds
  # starts or ends. The two new exons from such a split will have the same id,
  # causing these edges to appear.
  ft <- ft[ft[,1] != ft[,2],]
  ft[ft[,1] == ft[,2],]

  rEG <- graph::ftM2graphNEL(unique(ft))

  # Color exons
  cols <- RColorBrewer::brewer.pal(4, "Paired")
  interestcolor <- list(
    "excludedA" = cols[1],
    "includedA" = cols[2],
    "excludedB" = cols[3],
    "includedB" = cols[4])

  # Create a list of exon names with colors
  listA <- rep(interestcolor[2], length(exonnamesA))
  names(listA) <- exonnamesA
  listB <- rep(interestcolor[4], length(exonnamesB))
  names(listB) <- exonnamesB

  nAttrs <- list()
  eAttrs <- list()
  nAttrs$fillcolor <- c(
    "start" = "gray",
    listA,
    listB,
    "stop" = "gray")

  Rgraphviz::plot(
    rEG,
    attrs = list(
      node = list(
        shape = "box",
        fontsize = 42),
      edge = list(),
      graph = list(rankdir = rankdir)), # left-to-right: LR, top-to-bottom: TB
    nodeAttrs = nAttrs,
    edgeAttrs = eAttrs)

}

.validatePlotFusionTranscriptsGraphParams <- function(
  fusion,
  edb,
  whichTranscripts,
  rankdir
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()

  # Check parameters
  argument_checker <- .is.fusion.valid(argument_checker, fusion)
  argument_checker <- .is.edb.valid(argument_checker, edb, fusion)
  argument_checker <- .is.whichTranscripts.valid(
    argument_checker,
    whichTranscripts,
    fusion)

  # The "withinExon" option for whichTranscripts will not work for the fusion transcripts graph plot. Check if that was
  # the wanted option
  if (whichTranscripts[[1]] == "withinExon") {
    ArgumentCheck::addError(
      msg = paste0("The \"withinExon\" option for whichTranscripts will not ",
                   "work for the fusion transcripts graph plot, as only ",
                   "complete exons can be shown."),
      argcheck = argument_checker
    )
  }

  if (class(rankdir) != "character" || !rankdir %in% c("LR", "TB")) {
    ArgumentCheck::addError(
      msg = "'rankdir' must be one of \"LR\" or \"TB\"",
      argcheck = argument_checker
    )
  }

  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)
}
