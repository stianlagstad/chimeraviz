#' Plot a fusion event with transcripts, coverage and ideograms.
#'
#' This function creates a plot with information about transcripts, coverage,
#' location and more.
#'
#' plotFusion() will dispatch to either plotFusionSeparate() or
#' plotFusionTogether(). plotFusionSeparate() will plot the fusion gene
#' partners in separate graphs shown next to each other, while
#' plotFusionTogether() will plot the fusion gene partners in the same graph
#' with the same x-axis. plotFusion() will dispatch to plotFusionTogether() if
#' the fusion gene partners are on the same strand, same chromosome and are
#' close together (<=50,000 bp apart).
#'
#' @param fusion The Fusion object to plot.
#' @param edb The ensembldb object that will be used to fetch data.
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
#' @return Creates a fusion plot.
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
#' png(pngFilename, width = 1000, height = 750)
#' # Plot!
#' plotFusion(
#'   fusion = fusion,
#'   bamfile = bamfile5267,
#'   edb = edb,
#'   nonUCSC = TRUE)
#' # Close device
#' dev.off()
#'
#' @export
plotFusion <- function(
  fusion,
  edb = NULL,
  bamfile,
  whichTranscripts = "exonBoundary",
  ylim = c(0, 1000),
  nonUCSC = TRUE,
  reduceTranscripts = FALSE) {

  fusion <- .validateFusionPlotParams(fusion, edb, bamfile, whichTranscripts)

  # Decide whether or not to plot the fusion on the same x-axis
  if (fusion@geneA@chromosome == fusion@geneB@chromosome) {
    message("Fusion is intrachromosomal..")
    if (as.character(fusion@geneA@strand) != as.character(fusion@geneB@strand)) {
      message("..but on different strands. Plot separate!")
      plotFusionSeparate(
        fusion = fusion,
        edb = edb,
        bamfile = bamfile,
        whichTranscripts = whichTranscripts,
        ylim = ylim,
        nonUCSC = nonUCSC,
        reduceTranscripts = reduceTranscripts)
    } else if (abs(fusion@geneA@breakpoint - fusion@geneB@breakpoint) > 50000) {
      message("..but too far apart. Plot separate!")
      plotFusionSeparate(
        fusion = fusion,
        edb = edb,
        bamfile = bamfile,
        whichTranscripts = whichTranscripts,
        ylim = ylim,
        nonUCSC = nonUCSC,
        reduceTranscripts = reduceTranscripts)
    } else {
      message("..and close together..")
      # If both genes are on the minus strand, and the location of geneA is to
      # the right of geneB, then we still want to plot them separatly, because
      # we want geneA (the downstream fusion gene partner) to be on the left
      # side of the plot.
      if (fusion@geneA@strand == "-" &&
          fusion@geneB@strand == "-" &&
          min(start(unlist(fusion@geneA@transcripts))) <
          min(start(unlist(fusion@geneB@transcripts)))) {
        message("..but since both genes are on the minus strand the x-axis will be flipped in the plot. And since geneA will then be located to \" the right of\" geneB, we will still plot this separatly. Plot separate!")
        plotFusionSeparate(
          fusion = fusion,
          edb = edb,
          bamfile = bamfile,
          whichTranscripts = whichTranscripts,
          ylim = ylim,
          nonUCSC = nonUCSC,
          reduceTranscripts = reduceTranscripts)
      } else {
        message("Plot together!")
        plotFusionTogether(
          fusion = fusion,
          edb = edb,
          bamfile = bamfile,
          whichTranscripts = whichTranscripts,
          ylim = ylim,
          nonUCSC = nonUCSC,
          reduceTranscripts = reduceTranscripts)
      }
    }
  } else {
    message("Fusion is interchromosomal. Plot separate!")
    plotFusionSeparate(
      fusion = fusion,
      edb = edb,
      bamfile = bamfile,
      whichTranscripts = whichTranscripts,
      ylim = ylim,
      nonUCSC = nonUCSC,
      reduceTranscripts = reduceTranscripts)
  }
}

#' @name plotFusion
#'
#' @export
plotFusionSeparate <- function(
  fusion,
  edb,
  bamfile,
  whichTranscripts = "exonBoundary",
  ylim = c(0, 1000),
  nonUCSC = TRUE,
  reduceTranscripts = FALSE) {

  fusion <- .validateFusionPlotParams(fusion, edb, bamfile, whichTranscripts)

  # Create tracks

  # Select transcripts
  transcriptsA <- selectTranscript(fusion@geneA, whichTranscripts)
  transcriptsB <- selectTranscript(fusion@geneB, whichTranscripts)

  # Check that there's at least one transcript that has the fusion breakpoint
  # within the transcript.
  .checkThatBreakpointsAreWithinTranscripts(fusion, transcriptsA, transcriptsB)

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
    trAtrack <- Gviz::GeneRegionTrack(
      data.frame(reducedTranscriptsA),
      chromosome = fusion@geneA@chromosome)
    trBtrack <- Gviz::GeneRegionTrack(
      data.frame(reducedTranscriptsB),
      chromosome = fusion@geneB@chromosome)
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
    showTitle = FALSE, # hide left panel
    showId = TRUE, # show transcript id
    col = "transparent",
    fontcolor.group = "black",
    just.group = "below",
    min.distance = 0,
    fontsize.group = 25,
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
    add53 = TRUE,
    add35 = TRUE,
    col = "black", # line color
    fontcolor = "black",
    fontsize = 15
  )

  # Create ideogram tracks

  # Read cytoband information depending on genome version
  if (fusion@genomeVersion == "hg19") {
    cytobandFile <- system.file(
      "extdata",
      "UCSC.HG19.Human.CytoBandIdeogram.txt",
      package="chimeraviz")
  } else {
    cytobandFile <- system.file(
      "extdata",
      "UCSC.HG38.Human.CytoBandIdeogram.txt",
      package="chimeraviz")
  }
  cytoband <- utils::read.table(cytobandFile)
  # Set names to what Gviz expects
  names(cytoband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

  idTrackA <- Gviz::IdeogramTrack(
    genome = fusion@genomeVersion,
    chromosome = fusion@geneA@chromosome,
    bands = cytoband)
  idTrackB <- Gviz::IdeogramTrack(
    genome = fusion@genomeVersion,
    chromosome = fusion@geneB@chromosome,
    bands = cytoband)
  # Set display paramters
  ideogramDisplayParams <- list(
    showBandId = TRUE,
    showId = TRUE,
    fontcolor = "black",
    cex = 1.5)
  Gviz::displayPars(idTrackA) <- ideogramDisplayParams
  Gviz::displayPars(idTrackB) <- ideogramDisplayParams

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
      name = "RNA coverage",
      ylim = ylim,
      importFunction = importFunctionNonUCSC)
  } else {
    alTrack <- Gviz::AlignmentsTrack(
      bamfile,
      isPaired = TRUE,
      genome = fusion@genomeVersion,
      name = "RNA coverage",
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
    cex.axis = .6
  )

  # Display parameters for highlight tracks
  displayParsHighlightTracks <- list(
    fill = "lightgrey",
    col = "lightgrey",
    showTitle = FALSE
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

  grTrackHighlightA <- Gviz::HighlightTrack(
    trackList = trAtrack,
    start = highlightStartA,
    end = highlightEndA,
    chromosome = fusion@geneA@chromosome)
  Gviz::displayPars(grTrackHighlightA) <- displayParsHighlightTracks

  grTrackHighlightB <- Gviz::HighlightTrack(
    trackList = trBtrack,
    start = highlightStartB,
    end = highlightEndB,
    chromosome = fusion@geneB@chromosome)
  Gviz::displayPars(grTrackHighlightB) <- displayParsHighlightTracks

  # Highlight coverage tracks

  alTrackHighlightA <- Gviz::HighlightTrack(
    trackList = alTrack,
    start = highlightStartA,
    end = highlightEndA,
    chromosome = fusion@geneA@chromosome)
  Gviz::displayPars(alTrackHighlightA) <- displayParsHighlightTracks

  alTrackHighlightB <- Gviz::HighlightTrack(
    trackList = alTrack,
    start = highlightStartB,
    end = highlightEndB,
    chromosome = fusion@geneB@chromosome)
  Gviz::displayPars(alTrackHighlightB) <- displayParsHighlightTracks

  # Do some plotting

  # Create the grid we're gonna use for the plot
  nf <- grid::grid.layout(
    nrow = 1,
    ncol = 2,
    widths = grid::unit(c(1, 1), "null"))

  # Open a new page to plot it
  grid::grid.newpage()

  # Divide the page according to the grid we created
  grid::pushViewport(grid::viewport(layout = nf))

  # How tall should the transcripts row be? This depends on how many transcripts
  # we want to show. If reduceTranscripts=TRUE, then we only want to display
  # one transcript per gene.
  if (reduceTranscripts) {
    transcriptsRowHeight <- 1
  } else {
    # It depends on max(amountOfTranscriptsA, amountOfTranscriptsB)
    transcriptsRowHeight <- max(
      length(transcriptsA),
      length(transcriptsB))
    # Set a max value
    if (transcriptsRowHeight > 10) {
      transcriptsRowHeight <- 10
    }
  }

  # If the range for the plot is too short, then GenomeAxisTrack will show ticks too close to each other. This doesn't look nice. A way to solve it is to extend the plot range if the transcript is too short:
  plotRange <- max(end(trAtrack)) - min(start(trAtrack))
  if (plotRange < 20000) {
    offset <- 2500
  } else {
    offset <- 0
  }

  # Open row 1, column 1
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  # Plot transcriptsA and coverage with highlight
  res1 <- Gviz::plotTracks(
    collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
    c(idTrackA, grTrackHighlightA, alTrackHighlightA, axisTrack),
    from = min(start(trAtrack))-offset-10000, # 10k added so that we can see all transcript names
    to = max(end(trAtrack))+offset+10000, # 10k added so that we can see all transcript names
    sizes = c(1, transcriptsRowHeight, 2, .5),
    add = TRUE,
    background.title = "transparent",
    chromosome = rep(fusion@geneA@chromosome, 4),
    # Plot reverse if gene is at minus strand
    reverseStrand = geneAatMinusStrand)
  # Close row 1, column 1
  grid::popViewport(1)

  # If the range for the plot is too short, then GenomeAxisTrack will show ticks too close to each other. This doesn't look nice. A way to solve it is to extend the plot range if the transcript is too short:
  plotRange <- max(end(trBtrack)) - min(start(trBtrack))
  if (plotRange < 20000) {
    offset <- 2500
  } else {
    offset <- 0
  }

  # Open row 1, column 2
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  # Plot transcriptsB and coverage with highlight
  res2 <- Gviz::plotTracks(
    collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
    c(idTrackB, grTrackHighlightB, alTrackHighlightB, axisTrack),
    from = min(start(trBtrack))-offset-10000, # 10k added so that we can see all transcript names
    to = max(end(trBtrack))+offset+10000, # 10k added so that we can see all transcript names
    sizes = c(1, transcriptsRowHeight, 2, .5),
    add = TRUE,
    background.title = "transparent",
    chromosome = rep(fusion@geneB@chromosome, 4),
    # Plot reverse if gene is at minus strand
    reverseStrand = geneBatMinusStrand)
  # Close row 1, column 2
  grid::popViewport(1)

  # Draw bezier curve if we're plotting exon boundary transcripts

  if (length(fusion@geneA@transcripts[mcols(fusion@geneA@transcripts)$transcriptCategory == "exonBoundary"]) > 0 &
      length(fusion@geneB@transcripts[mcols(fusion@geneB@transcripts)$transcriptCategory == "exonBoundary"])) {

    # We need to figure out which exons in trAtrack and trBtrack matches the
    # fusion breakpoints fusion@geneA@breakpoint and fusion@geneB@breakpoint:
    if (geneAatMinusStrand) {
      exonBreakpointA <- mcols(trAtrack[GenomicRanges::start(trAtrack@range) == fusion@geneA@breakpoint]@range)$exon
    } else {
      exonBreakpointA <- mcols(trAtrack[GenomicRanges::end(trAtrack@range) == fusion@geneA@breakpoint]@range)$exon
    }
    if (geneBatMinusStrand) {
      exonBreakpointB <- mcols(trBtrack[GenomicRanges::end(trBtrack@range) == fusion@geneB@breakpoint]@range)$exon
    } else {
      exonBreakpointB <- mcols(trBtrack[GenomicRanges::start(trBtrack@range) == fusion@geneB@breakpoint]@range)$exon
    }

    # Get coordinates for the exons in the plot
    exonCoordinates1 <- Gviz::coords(res1$GeneRegionTrack)
    exonCoordinates2 <- Gviz::coords(res2$GeneRegionTrack)

    # Since we know which exon in geneA meets the breakpoint, its easy to fetch
    # the coordinates for the breakpoint exon in geneA:
    exonA <- exonCoordinates1[exonBreakpointA, , drop=FALSE]
    # Do the same for geneB:
    exonB <- exonCoordinates2[exonBreakpointB, , drop=FALSE]

    # The coordinates we're really after is the top right corner of exonA (or
    # top left corner if on the minus strand), and the top left (top right if on
    # minus strand) corner of exon B
    if (geneAatMinusStrand) {
      # top left corner:
      xPosGeneA <- min(exonA[, 1])
      yPosGeneA <- min(exonA[, 2])
    } else {
      # top right corner:
      xPosGeneA <- min(exonA[, 3])
      yPosGeneA <- min(exonA[, 2])
    }
    if (geneBatMinusStrand) {
      # top right corner:
      xPosGeneB <- min(exonB[, 3])
      yPosGeneB <- min(exonB[, 2])
    } else {
      # top left corner:
      xPosGeneB <- min(exonB[, 1])
      yPosGeneB <- min(exonB[, 2])
    }

    # This will keep the coordinates even though the dev.size change:
    grid::pushViewport(
      grid::viewport(
        xscale = c(0, grDevices::dev.size(units="px")[1]), # x goes from 0 to device width
        yscale = c(grDevices::dev.size(units="px")[2], 0))) # y goes from device height to 0

    # Calculate the bezier curve width. The following will have the effect of
    # setting the line width to 10 if there are 50 or more reads that support the
    # fusion. If there are between 1 and 50 reads supporting the fusion, then the
    # curve width will be scaled accordingly.
    supportingReads <- fusion@spanningReadsCount+fusion@splitReadsCount
    curveWidth <- .scaleListToInterval(
      c(supportingReads,
        1,
        20),
      1,
      5)
    curveWidth <- curveWidth[1]

    # How far above the transcript should the bezier curve control points be?
    bezierControlPointOffset <- 18

    # The different transcripts will be plotted at different heights. Choose the
    # y position closest to the top of the plot. That is the lowest y value of
    # the two, since the y axis is flipped
    yPos <- min(yPosGeneA, yPosGeneB)

    # Draw the bezier curve between the transcripts
    grid::grid.bezier(
      rep(c(xPosGeneA, xPosGeneB), each = 2), # x positions
      c(yPos,
        yPos - bezierControlPointOffset,
        yPos - bezierControlPointOffset,
        yPos), # y positions
      default.units = "native",
      gp = grid::gpar(
        col = "red",
        lwd = curveWidth))

    # Where should the number of supporting reads be plotted? Calculate position
    numberOffset <- 7
    numberPositionX <- xPosGeneA + (xPosGeneB - xPosGeneA) / 2
    numberPositionY <- yPos - bezierControlPointOffset - numberOffset

    grid::grid.text(
      supportingReads,
      x = numberPositionX/grDevices::dev.size(units="px")[1],
      y = 1 - numberPositionY/grDevices::dev.size(units="px")[2],
      vp = grid::viewport(
        xscale = c(0, grDevices::dev.size(units="px")[1]),
        yscale = c(grDevices::dev.size(units="px")[2], 0)),
      gp = grid::gpar(
        fontsize = 20))

  }
}

#' @name plotFusion
#'
#' @export
plotFusionTogether <- function(
  fusion,
  edb,
  bamfile,
  whichTranscripts = "exonBoundary",
  ylim = c(0, 1000),
  nonUCSC = TRUE,
  reduceTranscripts = FALSE) {

  fusion <- .validateFusionPlotParams(fusion, edb, bamfile, whichTranscripts)

  # Create tracks

  # Select transcripts
  transcriptsA <- selectTranscript(fusion@geneA, whichTranscripts)
  transcriptsB <- selectTranscript(fusion@geneB, whichTranscripts)

  # Check that there's at least one transcript that has the fusion breakpoint
  # within the transcript.
  .checkThatBreakpointsAreWithinTranscripts(fusion, transcriptsA, transcriptsB)

  # Set seqlevels to UCSC
  GenomeInfoDb::seqlevelsStyle(transcriptsA) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(transcriptsB) <- "UCSC"
  if (reduceTranscripts) {
    reducedTranscriptsA <- GenomicRanges::reduce(unlist(transcriptsA))
    reducedTranscriptsB <- GenomicRanges::reduce(unlist(transcriptsB))
    # Add transcript metadata column to make Gviz group correctly
    mcols(reducedTranscriptsA)$transcript <- fusion@geneA@name
    mcols(reducedTranscriptsB)$transcript <- fusion@geneB@name
    # Add symbol metadata column to make Gviz group correctly
    mcols(reducedTranscriptsA)$symbol <- fusion@geneA@name
    mcols(reducedTranscriptsB)$symbol <- fusion@geneB@name
    # Create the track
    trTrackBoth <- Gviz::GeneRegionTrack(
      rbind(data.frame(reducedTranscriptsA), data.frame(reducedTranscriptsB)),
      chromosome = fusion@geneA@chromosome)
  } else {
    # Create the track
    trTrackBoth <- Gviz::GeneRegionTrack(
      rbind(data.frame(transcriptsA), data.frame(transcriptsB)),
      chromosome = fusion@geneA@chromosome)
  }
  # Set display parameters
  transcriptDisplayParams <- list(
    showTitle = FALSE, # hide left panel
    showId = TRUE, # show transcript id
    col = "transparent",
    fontcolor.group = "black",
    just.group = "below",
    min.distance = 0,
    fontsize.group = 25,
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
  Gviz::displayPars(trTrackBoth) <- transcriptDisplayParams

  # Create axis track
  axisTrack <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axisTrack) <- list(
    add53 = TRUE,
    add35 = TRUE,
    col = "black", # line color
    fontcolor = "black",
    fontsize = 15
  )

  # Create ideogram tracks

  # Read cytoband information depending on genome version
  if (fusion@genomeVersion == "hg19") {
    cytobandFile <- system.file(
      "extdata",
      "UCSC.HG19.Human.CytoBandIdeogram.txt",
      package="chimeraviz")
  } else {
    cytobandFile <- system.file(
      "extdata",
      "UCSC.HG38.Human.CytoBandIdeogram.txt",
      package="chimeraviz")
  }
  cytoband <- utils::read.table(cytobandFile)
  # Set names to what Gviz expects
  names(cytoband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

  idTrack <- Gviz::IdeogramTrack(
    genome = fusion@genomeVersion,
    chromosome = fusion@geneA@chromosome,
    bands = cytoband)
  # Set display paramters
  ideogramDisplayParams <- list(
    showBandId = TRUE,
    showId = TRUE,
    fontcolor = "black",
    cex = 1.5)
  Gviz::displayPars(idTrack) <- ideogramDisplayParams

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
      name = "RNA coverage",
      ylim = ylim,
      importFunction = importFunctionNonUCSC)
  } else {
    alTrack <- Gviz::AlignmentsTrack(
      bamfile,
      isPaired = TRUE,
      genome = fusion@genomeVersion,
      name = "RNA coverage",
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
    cex.axis = .6
  )

  # Display parameters for highlight tracks
  displayParsHighlightTracks <- list(
    fill = "lightgrey",
    col = "lightgrey",
    showTitle = FALSE
  )

  # Helper variables that decide which way to plot the transcript later
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
  if (reduceTranscripts) {
    geneAindexes <- Gviz::symbol(trTrackBoth) == fusion@geneA@name
    geneBindexes <- Gviz::symbol(trTrackBoth) == fusion@geneB@name
  } else {
    geneAindexes <- Gviz::symbol(trTrackBoth) %in% names(transcriptsA)
    geneBindexes <- Gviz::symbol(trTrackBoth) %in% names(transcriptsB)
  }

  # Set all colors to exluded
  Gviz::feature(trTrackBoth)[Gviz::feature(trTrackBoth) == "5utr" & geneAindexes] <- "5utr_excludedA"
  Gviz::feature(trTrackBoth)[Gviz::feature(trTrackBoth) == "protein_coding" & geneAindexes] <- "protein_coding_excludedA"
  Gviz::feature(trTrackBoth)[Gviz::feature(trTrackBoth) == "3utr" & geneAindexes] <- "3utr_excludedA"

  Gviz::feature(trTrackBoth)[Gviz::feature(trTrackBoth) == "5utr" & geneBindexes] <- "5utr_excludedB"
  Gviz::feature(trTrackBoth)[Gviz::feature(trTrackBoth) == "protein_coding" & geneBindexes] <- "protein_coding_excludedB"
  Gviz::feature(trTrackBoth)[Gviz::feature(trTrackBoth) == "3utr" & geneBindexes] <- "3utr_excludedB"

  if (geneAatMinusStrand) {
    message(paste("As ",
                  fusion@geneA@name,
                  " is on the minus strand, the plot for this gene will be ",
                  "reversed",
                  sep = ""))
    highlightStartA <- fusion@geneA@breakpoint
    highlightEndA <- max(start(trTrackBoth[geneAindexes]) + width(trTrackBoth[geneAindexes]) - 1)

    # Color the exons included/not included in the fusion
    Gviz::feature(trTrackBoth)[start(trTrackBoth) <= highlightEndA &
                                 end(trTrackBoth) >= highlightStartA &
                                 Gviz::feature(trTrackBoth) == "5utr_excludedA" &
                                 geneAindexes] <- "5utr_includedA"

    Gviz::feature(trTrackBoth)[start(trTrackBoth) <= highlightEndA &
                                 end(trTrackBoth) >= highlightStartA &
                                 Gviz::feature(trTrackBoth) == "protein_coding_excludedA" &
                                 geneAindexes] <- "protein_coding_includedA"

    Gviz::feature(trTrackBoth)[start(trTrackBoth) <= highlightEndA &
                                 end(trTrackBoth) >= highlightStartA &
                                 Gviz::feature(trTrackBoth) == "3utr_excludedA" &
                                 geneAindexes] <- "3utr_includedA"
  } else {
    highlightStartA <- min(start(trTrackBoth[geneAindexes]))
    highlightEndA <- fusion@geneA@breakpoint

    # Color the exons included/not included in the fusion
    Gviz::feature(trTrackBoth)[start(trTrackBoth) >= highlightStartA &
                                 end(trTrackBoth) <= highlightEndA &
                                 Gviz::feature(trTrackBoth) == "5utr_excludedA" &
                                 geneAindexes] <- "5utr_includedA"

    Gviz::feature(trTrackBoth)[start(trTrackBoth) >= highlightStartA &
                                 end(trTrackBoth) <= highlightEndA &
                                 Gviz::feature(trTrackBoth) == "protein_coding_excludedA" &
                                 geneAindexes] <- "protein_coding_includedA"

    Gviz::feature(trTrackBoth)[start(trTrackBoth) >= highlightStartA &
                                 end(trTrackBoth) <= highlightEndA &
                                 Gviz::feature(trTrackBoth) == "3utr_excludedA" &
                                 geneAindexes] <- "3utr_includedA"
  }

  if (geneBatMinusStrand) {
    message(paste("As ",
                  fusion@geneB@name,
                  " is on the minus strand, the plot for this gene will be ",
                  "reversed",
                  sep = ""))
    highlightStartB <- min(start(trTrackBoth[geneBindexes]))
    highlightEndB <- fusion@geneB@breakpoint

    # Color the exons included/not included in the fusion
    Gviz::feature(trTrackBoth)[start(trTrackBoth) <= highlightEndB &
                                 end(trTrackBoth) >= highlightStartB &
                                 Gviz::feature(trTrackBoth) == "5utr_excludedB" &
                                 geneBindexes] <- "5utr_includedB"

    Gviz::feature(trTrackBoth)[start(trTrackBoth) <= highlightEndB &
                                 end(trTrackBoth) >= highlightStartB &
                                 Gviz::feature(trTrackBoth) == "protein_coding_excludedB" &
                                 geneBindexes] <- "protein_coding_includedB"

    Gviz::feature(trTrackBoth)[start(trTrackBoth) <= highlightEndB &
                                 end(trTrackBoth) >= highlightStartB &
                                 Gviz::feature(trTrackBoth) == "3utr_excludedB" &
                                 geneBindexes] <- "3utr_includedB"
  } else {
    highlightStartB <- fusion@geneB@breakpoint
    highlightEndB <- max(start(trTrackBoth[geneBindexes]) + width(trTrackBoth[geneBindexes]) - 1)

    # Color the exons included/not included in the fusion
    Gviz::feature(trTrackBoth)[start(trTrackBoth) >= highlightStartB &
                                 end(trTrackBoth) <= highlightEndB &
                                 Gviz::feature(trTrackBoth) == "5utr_excludedB" &
                                 geneBindexes] <- "5utr_includedB"

    Gviz::feature(trTrackBoth)[start(trTrackBoth) >= highlightStartB &
                                 end(trTrackBoth) <= highlightEndB &
                                 Gviz::feature(trTrackBoth) == "protein_coding_excludedB" &
                                 geneBindexes] <- "protein_coding_includedB"

    Gviz::feature(trTrackBoth)[start(trTrackBoth) >= highlightStartB &
                                 end(trTrackBoth) <= highlightEndB &
                                 Gviz::feature(trTrackBoth) == "3utr_excludedB" &
                                 geneBindexes] <- "3utr_includedB"
  }

  # Update transcript track with colors
  Gviz::displayPars(trTrackBoth) <- interestcolor

  # Highlight transcript track

  grTrackHighlight <- Gviz::HighlightTrack(
    trackList = trTrackBoth,
    start = c(highlightStartA, highlightStartB),
    end = c(highlightEndA, highlightEndB),
    chromosome = fusion@geneA@chromosome)
  Gviz::displayPars(grTrackHighlight) <- displayParsHighlightTracks

  # Highlight coverage track

  alTrackHighlight <- Gviz::HighlightTrack(
    trackList = alTrack,
    start = c(highlightStartA, highlightStartB),
    end = c(highlightEndA, highlightEndB),
    chromosome = fusion@geneA@chromosome)
  Gviz::displayPars(alTrackHighlight) <- displayParsHighlightTracks

  # Do some plotting

  # Open a new page to plot
  grid::grid.newpage()

  # How tall should the transcripts row be? This depends on how many transcripts
  # we want to show. If reduceTranscripts=TRUE, then we only want to display
  # one transcript per gene.
  if (reduceTranscripts) {
    transcriptsRowHeight <- 1
  } else {
    # It depends on max(amountOfTranscriptsA, amountOfTranscriptsB)
    transcriptsRowHeight <- max(
      length(transcriptsA),
      length(transcriptsB))
    # Set a max value
    if (transcriptsRowHeight > 10) {
      transcriptsRowHeight <- 10
    }
  }

  # If the range for the plot is too short, then GenomeAxisTrack will show ticks too close to each other. This doesn't look nice. A way to solve it is to extend the plot range if the transcript is too short:
  plotRange <- max(end(trTrackBoth)) - min(start(trTrackBoth))
  if (plotRange < 20000) {
    offset <- 2500
  } else {
    offset <- 0
  }

  # Plot transcripts and coverage with highlight
  res <- Gviz::plotTracks(
    collapse = FALSE, # without this gviz create cluster_X entries in the GeneRegionTrack
    c(idTrack, grTrackHighlight, alTrackHighlight, axisTrack),
    from = min(start(trTrackBoth))-offset-10000, # 10k added so that we can see all transcript names
    to = max(end(trTrackBoth))+offset+10000, # 10k added so that we can see all transcript names
    sizes = c(1, transcriptsRowHeight, 2, 1),
    add = TRUE,
    background.title = "transparent",
    chromosome = rep(fusion@geneA@chromosome, 4),
    # Plot reverse if gene is at minus strand
    reverseStrand = geneAatMinusStrand)

  # Draw bezier curve if we're plotting exon boundary transcripts

  if (length(fusion@geneA@transcripts[mcols(fusion@geneA@transcripts)$transcriptCategory == "exonBoundary"]) > 0 &
      length(fusion@geneB@transcripts[mcols(fusion@geneB@transcripts)$transcriptCategory == "exonBoundary"]) > 0) {

    # We need to figure out which exons in trTrackBoth matches the
    # fusion breakpoints fusion@geneA@breakpoint and fusion@geneB@breakpoint:
    if (geneAatMinusStrand) {
      exonBreakpointA <- mcols(trTrackBoth[GenomicRanges::start(trTrackBoth@range) == fusion@geneA@breakpoint]@range)$exon
    } else {
      exonBreakpointA <- mcols(trTrackBoth[GenomicRanges::end(trTrackBoth@range) == fusion@geneA@breakpoint]@range)$exon
    }
    if (geneBatMinusStrand) {
      exonBreakpointB <- mcols(trTrackBoth[GenomicRanges::end(trTrackBoth@range) == fusion@geneB@breakpoint]@range)$exon
    } else {
      exonBreakpointB <- mcols(trTrackBoth[GenomicRanges::start(trTrackBoth@range) == fusion@geneB@breakpoint]@range)$exon
    }

    # Get coordinates for the exons in the plot
    exonCoordinates <- Gviz::coords(res$GeneRegionTrack)

    # Since we know which exon in geneA meets the breakpoint, its easy to fetch
    # the coordinates for the breakpoint exon in geneA:
    exonA <- exonCoordinates[exonBreakpointA, , drop=FALSE]
    # Do the same for geneB:
    exonB <- exonCoordinates[exonBreakpointB, , drop=FALSE]

    # The coordinates we're really after is the top right corner of exonA (or
    # top left corner if on the minus strand), and the top left (top right if on
    # minus strand) corner of exon B
    if (geneAatMinusStrand) {
      # top left corner:
      xPosGeneA <- min(exonA[, 1])
      yPosGeneA <- min(exonA[, 2])
    } else {
      # top right corner:
      xPosGeneA <- min(exonA[, 3])
      yPosGeneA <- min(exonA[, 2])
    }
    if (geneBatMinusStrand) {
      # top right corner:
      xPosGeneB <- min(exonB[, 3])
      yPosGeneB <- min(exonB[, 2])
    } else {
      # top left corner:
      xPosGeneB <- min(exonB[, 1])
      yPosGeneB <- min(exonB[, 2])
    }

    # This will keep the coordinates even though the dev.size change:
    grid::pushViewport(
      grid::viewport(
        xscale = c(0, grDevices::dev.size(units="px")[1]), # x goes from 0 to device width
        yscale = c(grDevices::dev.size(units="px")[2], 0))) # y goes from device height to 0

    # # Draw some point just to check that we have the right coordinates:
    # grid.points(topCornerExonA[1], topCornerExonA[2])
    # grid.points(topCornerExonB[1], topCornerExonB[2])

    # Calculate the bezier curve width. The following will have the effect of
    # setting the line width to 10 if there are 50 or more reads that support the
    # fusion. If there are between 1 and 50 reads supporting the fusion, then the
    # curve width will be scaled accordingly.
    supportingReads <- fusion@spanningReadsCount+fusion@splitReadsCount
    curveWidth <- .scaleListToInterval(
      c(supportingReads,
        1,
        20),
      1,
      5)
    curveWidth <- curveWidth[1]

    # How far above the transcript should the bezier curve control points be?
    bezierControlPointOffset <- 18

    # The different transcripts will be plotted at different heights. Choose the
    # y position closest to the top of the plot. That is the lowest y value of
    # the two, since the y axis is flipped
    yPos <- min(yPosGeneA, yPosGeneB)

    # Draw the bezier curve between the transcripts
    grid::grid.bezier(
      rep(c(xPosGeneA, xPosGeneB), each = 2), # x positions
      c(yPos,
        yPos - bezierControlPointOffset,
        yPos - bezierControlPointOffset,
        yPos), # y positions
      default.units = "native",
      gp = grid::gpar(
        col = "red",
        lwd = curveWidth))

    # Where should the number of supporting reads be plotted? Calculate position
    numberOffset <- 7
    numberPositionX <- xPosGeneA + (xPosGeneB - xPosGeneA) / 2
    numberPositionY <- yPos - bezierControlPointOffset - numberOffset

    grid::grid.text(
      supportingReads,
      x = numberPositionX/grDevices::dev.size(units="px")[1],
      y = 1 - numberPositionY/grDevices::dev.size(units="px")[2],
      vp = grid::viewport(
        xscale = c(0, grDevices::dev.size(units="px")[1]),
        yscale = c(grDevices::dev.size(units="px")[2], 0)),
      gp = grid::gpar(
        fontsize = 20))

  }
}

#' Alternative import function for Gviz::AlignmentsTrack
#'
#' This alternative import function for use with Gviz::AlignmentsTrack imports
#' a bamfile with non-UCSC chromosome names.
#'
#' @param file The bamfile.
#' @param selection Which regions to get from the bamfile.
#'
#' @return A GRanges object with coverage data for the selection.
importFunctionNonUCSC <- function (file, selection) {

  indNames <- c(sub("\\.bam$", ".bai", file), paste(file,
                                                    "bai", sep = "."))
  index <- NULL
  for (i in indNames) {
    if (file.exists(i)) {
      index <- i
      break
    }
  }
  if (is.null(index))
    stop("Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t",
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  pairedEnd <- parent.env(environment())[["._isPaired"]]
  if (is.null(pairedEnd))
    pairedEnd <- TRUE
  bf <- Rsamtools::BamFile(file, index = index, asMates = pairedEnd)
  # Rename chromosome names from UCSD to non-UCSC
  GenomeInfoDb::seqlevels(selection) <- sub(
    "chr",
    "",
    GenomeInfoDb::seqlevels(selection))
  param <- Rsamtools::ScanBamParam(
    which = selection,
    what = Rsamtools::scanBamWhat(),
    tag = "MD",
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
  reads <- if (as.character(seqnames(selection)[1]) %in% names(Rsamtools::scanBamHeader(bf)$targets))
    Rsamtools::scanBam(bf, param = param)[[1]]
  else list()
  md <- if (is.null(reads$tag$MD))
    rep(as.character(NA), length(reads$pos))
  else reads$tag$MD
  if (length(reads$pos)) {
    layed_seq <- GenomicAlignments::sequenceLayer(reads$seq, reads$cigar)
    region <- unlist(Rsamtools::bamWhich(param), use.names = FALSE)
    ans <- stackStrings(
      layed_seq,
      start(region),
      end(region),
      shift = reads$pos - 1L,
      Lpadding.letter = "+",
      Rpadding.letter = "+")
    names(ans) <- seq_along(reads$qname)
  }
  else {
    ans <- DNAStringSet()
  }
  return(GRanges(
    # Return chromosome name prefixed with "chr"
    seqnames = if (is.null(reads$rname)) character() else paste("chr", reads$rname, sep=""),
    strand = if (is.null(reads$strand)) character() else reads$strand,
    ranges = IRanges(start = reads$pos, width = reads$qwidth),
    id = reads$qname, cigar = reads$cigar, mapq = reads$mapq,
    flag = reads$flag, md = md, seq = ans, isize = reads$isize,
    groupid = if (pairedEnd) reads$groupid else seq_along(reads$pos),
    status = if (pairedEnd) reads$mate_status else rep(factor("unmated", levels = c("mated", "ambiguous", "unmated")), length(reads$pos))))
}
