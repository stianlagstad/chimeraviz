#' Plot possible fusion transcripts based on annotation.
#'
#' This function takes a fusion object and an ensembldb object and plots the
#' reduced version of the fusion transcript. This transcript consist of the
#' "mashed together" version of all possible fusion transcripts based on known
#' annotations. See plotFusionTranscripts() for the not-reduced version of this
#' plot. If a bamfile is specified, the fusion transcript will be plotted with
#' coverage information.
#'
#' Note that the transcript database used (the edb object) must have the same
#' seqnames as any bamfile used. Otherwise the coverage data will be wrong.
#'
#' @param fusion The Fusion object to plot.
#' @param edb The edb object that will be used to fetch data.
#' @param bamfile The bamfile with RNA-seq data.
#' @param whichTranscripts This character vector decides which transcripts are
#' to be plotted. Can be "exonBoundary", "withinExon", "withinIntron",
#' "intergenic", or a character vector with specific transcript ids. Default
#' value is "exonBoundary".
#'
#' @return Creates a fusion transcript plot.
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
#' # bamfile with reads in the regions of this fusion event
#' bamfile5267 <- system.file(
#'   "extdata",
#'   "fusion5267and11759reads.bam",
#'   package="chimeraviz")
#' # Temporary file to store the plot
#' pngFilename <- tempfile(
#'   pattern = "fusionPlot",
#'   fileext = ".png",
#'   tmpdir = tempdir())
#' # Open device
#' png(pngFilename, width = 500, height = 500)
#' # Plot!
#' plotFusionTranscript(
#'   fusion = fusion,
#'   bamfile = bamfile5267,
#'   edb = edb)
#' # Close device
#' dev.off()
#'
#' @export
plotFusionTranscript <- function(
  fusion,
  edb = NULL,
  bamfile = NULL,
  whichTranscripts = "exonBoundary") {

  fusion <- .validateFusionPlotParams(fusion, edb, bamfile, whichTranscripts)

  # Select which transcripts to use
  transcriptsA <- selectTranscript(fusion@geneA, whichTranscripts)
  transcriptsB <- selectTranscript(fusion@geneB, whichTranscripts)

  # Check that there's at least one transcript that has the fusion breakpoint
  # within the transcript.
  .checkThatBreakpointsAreWithinTranscripts(fusion, transcriptsA, transcriptsB)

  # Reduce transcript into one for each gene
  transcriptsA <- reduce(unlist(transcriptsA))
  transcriptsB <- reduce(unlist(transcriptsB))
  # Add transcript metadata column to make Gviz group correctly
  mcols(transcriptsA)$transcript <- fusion@geneA@name
  mcols(transcriptsB)$transcript <- fusion@geneB@name
  # Add symbol metadata column to make Gviz show ids
  mcols(transcriptsA)$symbol <- fusion@geneA@name
  mcols(transcriptsB)$symbol <- fusion@geneB@name

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

  # Reverse order of exons if on minus strand
  if (fusion@geneA@strand == "-") {
    fusionExonsA <- rev(fusionExonsA)
  }
  if (fusion@geneB@strand == "-") {
    fusionExonsB <- rev(fusionExonsB)
  }

  # Build fusion transcript
  fusionTranscript <- append(
    fusionExonsA,
    fusionExonsB
  )

  # If we've got a bamfile, calculate coverage
  if (!is.null(bamfile)) {
    # Get coverage
    cov <- coverage(bamfile, param = Rsamtools::ScanBamParam(which = fusionTranscript))
    # Coverage only for my transcript
    cov <- cov[fusionTranscript]
    # Unlist
    cov <- suppressWarnings(unlist(cov))
    # Turn into numeric vector
    cov <- as.numeric(cov)

    # Create coverage track
    dTrack <- DataTrack(
      start = 1:length(cov),
      width = 1,
      chromosome = "chrNA",
      genome = "hg19",
      name = "Coverage",
      type = "h",
      col = 'orange',
      fill = 'orange',
      data = cov)
    # Set display parameters
    Gviz::displayPars(dTrack) <- list(
      showTitle = FALSE, # hide name of track
      background.panel = "transparent", # background color of the content panel
      background.title = "transparent", # background color for the title panels
      cex.axis = .6,
      cex.title = .6,
      col.axis = "black",
      col.coverage = "black",
      col.title = "black",
      coverageHeight = 0.08,
      fill.coverage = "orange",
      fontsize = 15,
      lty = 1,
      lty.coverage = 1,
      lwd = 0.5
    )
  }

  # Move the first exon down to position 1, and then remove the spaces between
  # all the others. Effecticely pushing all exons to the leftmost edge.
  for (i in 1:length(fusionTranscript)) {
    if (i == 1) {
      fusionTranscript[i] <- IRanges::shift(fusionTranscript[i], shift = 1-start(fusionTranscript[i]))
    } else {
      shiftDownBy <- -(start(fusionTranscript[i])-end(fusionTranscript[i-1]))+1
      fusionTranscript[i] <- IRanges::shift(fusionTranscript[i], shift = shiftDownBy)
    }
  }

  # Create new GRanges so that we can set the right seqname
  gr <- GRanges(seqnames = "chrNA", ranges = ranges(fusionTranscript))
  mcols(gr)$symbol <- mcols(fusionTranscript)$symbol

  # Create transcript track
  trTrack <- Gviz::GeneRegionTrack(
    gr,
    chromosome = "chrNA")
  # Set display parameters
  Gviz::displayPars(trTrack) <- list(
    background.panel = "transparent", # background color of the content panel
    background.title = "transparent", # background color for the title panels
    col.title = "black",
    showTitle = TRUE, # hide left panel
    showId = FALSE, # show transcript id
    col = "lightgray", # border around exons
    showTitle = FALSE # hide title
  )

  # Color exons
  cols <- RColorBrewer::brewer.pal(4, "Paired")
  interestcolor <- list(
    "excludedA" = cols[1],
    "includedA" = cols[2],
    "excludedB" = cols[3],
    "includedB" = cols[4])
  Gviz::feature(trTrack)[Gviz::symbol(trTrack) == fusion@geneA@name] <- "includedA"
  Gviz::feature(trTrack)[Gviz::symbol(trTrack) == fusion@geneB@name] <- "includedB"
  Gviz::displayPars(trTrack) <- interestcolor

  # Create axis track
  axisTrack <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axisTrack) <- list(
    col = "black", # line color
    fontcolor = "black",
    fontsize= 8,
    lwd = 1.5
  )

  # If we've got a bamfile, plot with coverage
  if (!is.null(bamfile)) {
    Gviz::plotTracks(
      list(dTrack, trTrack, axisTrack),
      main = paste(fusion@geneA@name, fusion@geneB@name, sep = ":"),
      chromosome = "chrNA")
  } else {
    Gviz::plotTracks(
      list(trTrack, axisTrack),
      main = paste(fusion@geneA@name, fusion@geneB@name, sep = ":"),
      chromosome = "chrNA")
  }

}
