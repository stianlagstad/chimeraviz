#' Plot transcripts for each partner gene in a fusion event.
#'
#' This function takes a fusion object and an ensembldb object and plots
#' transcripts for both genes, showing which parts of each genes are included
#' in the fusion event. If the bamfile parameter is set, then the coverage is
#' plotted beneath the transcripts.
#'
#' @param fusion The Fusion object to plot.
#' @param edb The edb object that will be used to fetch data.
#' @param bamfile The bamfile with RNA-seq data.
#' @param whichTranscripts This character vector decides which transcripts are
#' to be plotted. Can be "exonBoundary", "withinExon", "withinIntron",
#' "intergenic", or a character vector with specific transcript ids. Default
#' value is "exonBoundary".
#' @param ylim Limits for the coverage y-axis.
#' @param nonUCSC Boolean indicating whether or not the bamfile used has UCSC-
#' styled chromosome names (i.e. with the "chr" prefix). Setting this to true
#' lets you use a bamfile with chromosome names like "1" and "X", instead of
#' "chr1" and "chrX".
#' @param reduceTranscripts Boolean indicating whether or not to reduce all
#' transcripts into a single transcript for each partner gene.
#'
#' @return Creates a fusion transcripts plot.
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
#' plotTranscripts(
#'   fusion = fusion,
#'   edb = edb,
#'   bamfile = bamfile5267,
#'   nonUCSC = TRUE)
#' # Close device
#' dev.off()
#'
#' @importFrom ensembldb EnsDb
#'
#' @export
plotTranscripts <- function(
  fusion,
  edb = NULL,
  bamfile = NULL,
  whichTranscripts = "exonBoundary",
  nonUCSC = TRUE,
  ylim = c(0,1000),
  reduceTranscripts = FALSE) {

  .validatePlotFusionTranscriptsParams(fusion, edb, bamfile, whichTranscripts, nonUCSC, ylim, reduceTranscripts)
  fusion <- .getTranscriptsIfNotThere(fusion, edb)

  # Select which transcripts to use
  transcriptsA <- selectTranscript(fusion@geneA, whichTranscripts)
  transcriptsB <- selectTranscript(fusion@geneB, whichTranscripts)

  # Set seqlevels to UCSC
  GenomeInfoDb::seqlevelsStyle(transcriptsA) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(transcriptsB) <- "UCSC"

  if (reduceTranscripts) {
    reducedTranscriptsA <- GenomicRanges::reduce(unlist(transcriptsA))
    reducedTranscriptsB <- GenomicRanges::reduce(unlist(transcriptsB))
    # Add transcript metadata column to make Gviz group correctly
    mcols(reducedTranscriptsA)$transcript <- fusion@geneA@name
    mcols(reducedTranscriptsB)$transcript <- fusion@geneB@name
    # Create the tracks
    trAtrack <- Gviz::GeneRegionTrack(reducedTranscriptsA, chromosome = fusion@geneA@chromosome)
    trBtrack <- Gviz::GeneRegionTrack(reducedTranscriptsB, chromosome = fusion@geneB@chromosome)
    # Set symbol names on the tracks
    symbol(trAtrack) <- fusion@geneA@name
    symbol(trBtrack) <- fusion@geneB@name
  } else {
    # Create the tracks
    trAtrack <- Gviz::GeneRegionTrack(
      data.frame(transcriptsA),
      chromosome = fusion@geneA@chromosome)
    trBtrack <- Gviz::GeneRegionTrack(
      data.frame(transcriptsB),
      chromosome = fusion@geneB@chromosome)
  }

  # Set display parameters
  transcriptDisplayParams <- list(
    background.title = "transparent", # hide left panel
    showId = TRUE, # show transcript id
    col = "transparent",
    cex.group = 0.8,
    fontcolor.group = "black",
    just.group = "below",
    min.distance = 0,
    thinBoxFeature = c(
      "5utr_excludedA",
      "5utr_includedA",
      "3utr_excludedA",
      "3utr_includedA",

      "5utr_excludedB",
      "5utr_includedB",
      "3utr_excludedB",
      "3utr_includedB"
    ))
  Gviz::displayPars(trAtrack) <- transcriptDisplayParams
  Gviz::displayPars(trBtrack) <- transcriptDisplayParams

  # Create axis track
  axisTrack <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axisTrack) <- list(
    labelPos = "below",
    cex = 1,
    col = "black", # line color
    fontcolor = "black",
    lwd = 1.5
  )

  # Helper variables
  geneAatMinusStrand <- fusion@geneA@strand == "-"
  geneBatMinusStrand <- fusion@geneB@strand == "-"

  # Color options for exons included/excluded in fusion
  cols <- RColorBrewer::brewer.pal(4, "Paired")
  interestcolor <- list(
    "5utr_excludedA" = cols[1],
    "protein_coding_excludedA" = cols[1],
    "3utr_excludedA" = cols[1],

    "5utr_includedA" = cols[2],
    "protein_coding_includedA" = cols[2],
    "3utr_includedA" = cols[2],

    "5utr_excludedB" = cols[3],
    "protein_coding_excludedB" = cols[3],
    "3utr_excludedB" = cols[3],

    "5utr_includedB" = cols[4],
    "protein_coding_includedB" = cols[4],
    "3utr_includedB" = cols[4])

  # Set all colors to exluded
  if (reduceTranscripts) {
    Gviz::feature(trAtrack) <- "protein_coding_excludedA"
    Gviz::feature(trBtrack) <- "protein_coding_excludedB"
  } else {
    Gviz::feature(trAtrack)[Gviz::feature(trAtrack) == "5utr"] <- "5utr_excludedA"
    Gviz::feature(trAtrack)[Gviz::feature(trAtrack) == "protein_coding"] <- "protein_coding_excludedA"
    Gviz::feature(trAtrack)[Gviz::feature(trAtrack) == "3utr"] <- "3utr_excludedA"

    Gviz::feature(trBtrack)[Gviz::feature(trBtrack) == "5utr"] <- "5utr_excludedB"
    Gviz::feature(trBtrack)[Gviz::feature(trBtrack) == "protein_coding"] <- "protein_coding_excludedB"
    Gviz::feature(trBtrack)[Gviz::feature(trBtrack) == "3utr"] <- "3utr_excludedB"
  }

  if (geneAatMinusStrand) {
    message(paste("As ",
                  fusion@geneA@name,
                  " is on the minus strand, the plot for this gene will be ",
                  "reversed",
                  sep = ""))
    highlightStartA <- fusion@geneA@breakpoint
    highlightEndA <- max(start(trAtrack) + width(trAtrack) - 1)

    # Color the exons included/not included in the fusion
    Gviz::feature(trAtrack)[start(trAtrack) <= highlightEndA &
                              end(trAtrack) >= highlightStartA &
                              Gviz::feature(trAtrack) == "5utr_excludedA"] <- "5utr_includedA"

    Gviz::feature(trAtrack)[start(trAtrack) <= highlightEndA &
                              end(trAtrack) >= highlightStartA &
                              Gviz::feature(trAtrack) == "protein_coding_excludedA"] <- "protein_coding_includedA"

    Gviz::feature(trAtrack)[start(trAtrack) <= highlightEndA &
                              end(trAtrack) >= highlightStartA &
                              Gviz::feature(trAtrack) == "3utr_excludedA"] <- "3utr_includedA"

  } else {
    highlightStartA <- min(start(trAtrack))
    highlightEndA <- fusion@geneA@breakpoint

    # Color the exons included/not included in the fusion
    Gviz::feature(trAtrack)[start(trAtrack) >= highlightStartA &
                              end(trAtrack) <= highlightEndA &
                              Gviz::feature(trAtrack) == "5utr_excludedA"] <- "5utr_includedA"

    Gviz::feature(trAtrack)[start(trAtrack) >= highlightStartA &
                              end(trAtrack) <= highlightEndA &
                              Gviz::feature(trAtrack) == "protein_coding_excludedA"] <- "protein_coding_includedA"

    Gviz::feature(trAtrack)[start(trAtrack) >= highlightStartA &
                              end(trAtrack) <= highlightEndA &
                              Gviz::feature(trAtrack) == "3utr_excludedA"] <- "3utr_includedA"

  }

  if (geneBatMinusStrand) {
    message(paste("As ",
                  fusion@geneB@name,
                  " is on the minus strand, the plot for this gene will be ",
                  "reversed",
                  sep = ""))
    highlightStartB <- min(start(trBtrack))
    highlightEndB <- fusion@geneB@breakpoint

    # Color the exons included/not included in the fusion
    Gviz::feature(trBtrack)[start(trBtrack) <= highlightEndB &
                              end(trBtrack) >= highlightStartB &
                              Gviz::feature(trBtrack) == "5utr_excludedB"] <- "5utr_includedB"

    Gviz::feature(trBtrack)[start(trBtrack) <= highlightEndB &
                              end(trBtrack) >= highlightStartB &
                              Gviz::feature(trBtrack) == "protein_coding_excludedB"] <- "protein_coding_includedB"

    Gviz::feature(trBtrack)[start(trBtrack) <= highlightEndB &
                              end(trBtrack) >= highlightStartB &
                              Gviz::feature(trBtrack) == "3utr_excludedB"] <- "3utr_includedB"

  } else {
    highlightStartB <- fusion@geneB@breakpoint
    highlightEndB <- max(start(trBtrack) + width(trBtrack) - 1)

    # Color the exons included/not included in the fusion
    Gviz::feature(trBtrack)[start(trBtrack) >= highlightStartB &
                              end(trBtrack) <= highlightEndB &
                              Gviz::feature(trBtrack) == "5utr_excludedB"] <- "5utr_includedB"

    Gviz::feature(trBtrack)[start(trBtrack) >= highlightStartB &
                              end(trBtrack) <= highlightEndB &
                              Gviz::feature(trBtrack) == "protein_coding_excludedB"] <- "protein_coding_includedB"

    Gviz::feature(trBtrack)[start(trBtrack) >= highlightStartB &
                              end(trBtrack) <= highlightEndB &
                              Gviz::feature(trBtrack) == "3utr_excludedB"] <- "3utr_includedB"

  }

  # Update transcript track with colors
  Gviz::displayPars(trAtrack) <- interestcolor
  Gviz::displayPars(trBtrack) <- interestcolor

  # Highlight transcript tracks

  # Display parameters for highlight tracks
  displayParsHighlightTracks <- list(
    fill = "lightgrey",
    col = "lightgrey"
  )

  # Only plot highlight tracks if we actually have transcripts that has the
  # fusion breakpoint within them

  if (any(start(trAtrack) < fusion@geneA@breakpoint) &&
      any(fusion@geneA@breakpoint < end(trAtrack))) {

    grTrackHighlightA <- Gviz::HighlightTrack(
      trackList = trAtrack,
      start = highlightStartA,
      end = highlightEndA,
      chromosome = fusion@geneA@chromosome)
    Gviz::displayPars(grTrackHighlightA) <- displayParsHighlightTracks
  }

  if (any(start(trBtrack) < fusion@geneB@breakpoint) &&
      any(fusion@geneB@breakpoint < end(trBtrack))) {

    grTrackHighlightB <- Gviz::HighlightTrack(
      trackList = trBtrack,
      start = highlightStartB,
      end = highlightEndB,
      chromosome = fusion@geneB@chromosome)
    Gviz::displayPars(grTrackHighlightB) <- displayParsHighlightTracks
  }

  # Plot coverage?
  if (!is.null(bamfile)) {
    # Create alignment track
    if (nonUCSC) {
      # If the bam file has non-ucsc chromosome names, i.e. "1" instead of "chr1",
      # then we need to use a custom import function that adds "chr" to the
      # chromosome names, to keep Gviz happy. Gviz strongly prefers having the
      # "chr" prefix.
      alTrack <- Gviz::AlignmentsTrack(
        bamfile,
        isPaired = TRUE,
        genome = fusion@genomeVersion,
        ylim = ylim,
        importFunction = importFunctionNonUCSC)
    } else {
      alTrack <- Gviz::AlignmentsTrack(
        bamfile,
        isPaired = TRUE,
        genome = fusion@genomeVersion,
        ylim = ylim)
    }
    # Set display paramters
    Gviz::displayPars(alTrack) <- list(
      showTitle = FALSE, # hide name of track
      background.panel = "transparent", # background color of the content panel
      background.title = "transparent", # background color for the title panels
      col.axis = "black",
      col.coverage = "black",
      col.title = "black",
      fill.coverage = "orange",
      fontsize = 15,
      lwd.coverage = 0.4,
      type = "coverage",
      cex.title = .8,
      cex.axis = .6,
      show
    )

    # Highlight coverage tracks

    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(trAtrack) < fusion@geneA@breakpoint) &&
        any(fusion@geneA@breakpoint < end(trAtrack))) {

      alTrackHighlightA <- Gviz::HighlightTrack(
        trackList = alTrack,
        start = highlightStartA,
        end = highlightEndA,
        chromosome = fusion@geneA@chromosome)
      Gviz::displayPars(alTrackHighlightA) <- displayParsHighlightTracks
    }

    if (any(start(trBtrack) < fusion@geneB@breakpoint) &&
        any(fusion@geneB@breakpoint < end(trBtrack))) {

      alTrackHighlightB <- Gviz::HighlightTrack(
        trackList = alTrack,
        start = highlightStartB,
        end = highlightEndB,
        chromosome = fusion@geneB@chromosome)
      Gviz::displayPars(alTrackHighlightB) <- displayParsHighlightTracks
    }
  }

  # Do some plotting

  # Calculate heights in the plot grid according to how many transcripts there
  # are in each partnergene
  if (reduceTranscripts) {
    numberOfTranscriptsA <- 1
    numberOfTranscriptsB <- 1
  } else {
    numberOfTranscriptsA <- length(transcriptsA)
    numberOfTranscriptsB <- length(transcriptsB)
  }

  # Create the grid we're gonna use for the plot
  nf <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(numberOfTranscriptsA, numberOfTranscriptsB), "null"))

  # Open a new page to plot it
  grid::grid.newpage()

  # Divide the page according to the grid we created
  grid::pushViewport(grid::viewport(layout = nf))

  # Open row 1, column 1
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  # Plot transcriptsA and coverage with highlight
  # Plot coverage?
  if (!is.null(bamfile)) {

    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(trAtrack) < fusion@geneA@breakpoint) &&
        any(fusion@geneA@breakpoint < end(trAtrack))) {

      Gviz::plotTracks(
        collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
        list(grTrackHighlightA, alTrackHighlightA),
        sizes = c(5, 2),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = geneAatMinusStrand)
    } else {
      Gviz::plotTracks(
        collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
        list(trAtrack, alTrack),
        sizes = c(5, 2),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = geneAatMinusStrand)
    }
  } else {
    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(trAtrack) < fusion@geneA@breakpoint) &&
        any(fusion@geneA@breakpoint < end(trAtrack))) {

      Gviz::plotTracks(
        collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
        grTrackHighlightA,
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = geneAatMinusStrand)
    } else {
      Gviz::plotTracks(
        collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
        trAtrack,
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = geneAatMinusStrand)
    }
  }
  # Close row 1, column 1
  grid::popViewport(1)

  # Open row 1, column 1
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  # Plot transcriptsA and coverage with highlight
  # Plot coverage?
  if (!is.null(bamfile)) {
    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(trBtrack) < fusion@geneB@breakpoint) &&
        any(fusion@geneB@breakpoint < end(trBtrack))) {

      Gviz::plotTracks(
        collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
        list(grTrackHighlightB, alTrackHighlightB),
        sizes = c(5, 2),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = geneBatMinusStrand)
    } else {
      Gviz::plotTracks(
        collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
        list(trBtrack, alTrackHighlightB),
        sizes = c(5, 2),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = geneBatMinusStrand)
    }
  } else {
    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(trBtrack) < fusion@geneB@breakpoint) &&
        any(fusion@geneB@breakpoint < end(trBtrack))) {

      Gviz::plotTracks(
        collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
        grTrackHighlightB,
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = geneBatMinusStrand)
    } else {
      Gviz::plotTracks(
        collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
        trAtrack,
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = geneBatMinusStrand)
    }
  }
  # Close row 1, column 1
  grid::popViewport(1)
}

.validatePlotFusionTranscriptsParams <- function(
  fusion,
  edb,
  bamfile,
  whichTranscripts,
  nonUCSC,
  ylim,
  reduceTranscripts
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()

  # Check parameters
  argument_checker <- .is.fusion.valid(argument_checker, fusion)
  argument_checker <- .is.edb.valid(argument_checker, edb, fusion)
  if (!is.null(bamfile)) {
    argument_checker <- .is.bamfile.valid(argument_checker, bamfile)
  }
  argument_checker <- .is.whichTranscripts.valid(
    argument_checker,
    whichTranscripts,
    fusion)
  argument_checker <- .is.ylim.valid(argument_checker, ylim)
  argument_checker <- .is.parameter.boolean(
    argument_checker,
    nonUCSC,
    "nonUCSC")
  argument_checker <- .is.parameter.boolean(
    argument_checker,
    reduceTranscripts,
    "reduceTranscripts")

  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)
}
