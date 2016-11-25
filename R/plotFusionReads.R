#' Create a plot of the reads supporting the given fusion.
#'
#' This function takes a Fusion object and plots the reads supporting the fusion
#' on top of the fusion sequence (fusion@junctionSequence), provided that
#' addFusionReadsAlignment() has been run earlier in order to add fusion reads
#' alignment data to the fusion object.
#'
#' Note that the package used for plotting, Gviz, is strict on chromosome names.
#' If the plot produced doesn't show the reads, the problem might be solved by
#' naming the fusion sequence "chrNA".
#'
#' @param fusion The Fusion object to plot.
#' @param showAllNucleotides Boolean indicating whether or not to show all
#' nucleotides. If FALSE, then only nucleotideAmount amount of nucleotides will
#' be shown on each end of the fusion junction. If TRUE, then the whole fusion
#' junction sequence will be shown.
#' @param nucleotideAmount The number of nucleotides to show on each end of the
#' fusion junction sequence. Defaults to 10. Only applicable if
#' showAllNucleotides is set to TRUE.
#'
#' @return Creates a fusion reads plot.
#'
#' @seealso addFusionReadsAlignment
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' # Find the specific fusion we have aligned reads for
#' fusion <- getFusionById(fusions, 5267)
#' bamfile <- system.file(
#'   "extdata",
#'   "5267readsAligned.bam",
#'   package="chimeraviz")
#' # Add the bam file of aligned fusion reads to the fusion object
#' fusion <- addFusionReadsAlignment(fusion, bamfile)
#' # Temporary file to store the plot
#' pngFilename <- tempfile(
#'   pattern = "fusionPlot",
#'   fileext = ".png",
#'   tmpdir = tempdir())
#' # Calculate image size based on supporting reads and lenght of junction
#' # sequence.
#' imageWidth <- (nchar(fusion@geneA@junctionSequence) + nchar(fusion@geneB@junctionSequence)) * 15
#' imageHeight <- (fusion@splitReadsCount+fusion@spanningReadsCount) * 20
#' # Open device
#' png(pngFilename, width = imageWidth, height = imageHeight)
#' # Now we can plot
#' plotFusionReads(fusion)
#' # Close device
#' dev.off()
#'
#' @export
plotFusionReads <- function(fusion, showAllNucleotides = TRUE, nucleotideAmount = 10) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  # Check if the Gviz::AlignmentsTrack object in
  # fusion@fusionReadsAlignment actually holds data, and is not just the empty
  # object it is after loading the initial fusion data. In other words, make
  # sure addFusionReadsAlignment() has been run.
  if (class(fusion@fusionReadsAlignment) != "ReferenceAlignmentsTrack") {
    stop(paste("fusion@fusionReadsAlignment should hold a",
               "Gviz::ReferenceAlignmentsTrack object.",
               "See addFusionReadsAlignment()."))
  }

  # Create Biostrings::DNAStringSet of fusion sequence
  fusionSequence <- Biostrings::DNAStringSet(
    x = c(fusion@geneA@junctionSequence,
          fusion@geneB@junctionSequence))
  names(fusionSequence) <- "chrNA"
  # And also of a "zoomed-in" version only displaying some nucleotides on each
  # end of the fusion junction
  fusionSequenceShort <- Biostrings::DNAStringSet(
    x = c(
      substr(
        fusion@geneA@junctionSequence,
        nchar(fusion@geneA@junctionSequence)-nucleotideAmount+1,
        nchar(fusion@geneA@junctionSequence)),
      substr(
        fusion@geneB@junctionSequence,
        nchar(fusion@geneB@junctionSequence)-nucleotideAmount+1 -1,
        # Added -1 again on the line above because Gviz does not show the last
        # nucleotide. Ref: http://permalink.gmane.org/gmane.science.biology.informatics.conductor.devel/5726
        # Or if that link is down, google
        # "Gviz doesn't provide enough horizontial space"
        nchar(fusion@geneB@junctionSequence))))
  names(fusionSequenceShort) <- "chrNA"

  # Load fusion sequences to tracks
  ref_chimera_short <- Gviz::SequenceTrack(
    fusionSequenceShort,
    genome = fusion@genomeVersion,
    name = "Fusion Sequence")
  ref_chimera <- Gviz::SequenceTrack(
    fusionSequence,
    genome = fusion@genomeVersion,
    name = "Fusion Sequence")

  # DISPLAY PARAMETERS

  # Set display parameters for the fusion sequence
  Gviz::displayPars(ref_chimera) <- list(
    showTitle = FALSE,
    col = "blue", # The color of the line when no indiviual letters can be plotted due to size limitations.
    fontsize = 10
  )
  Gviz::displayPars(ref_chimera_short) <- list(
    showTitle = FALSE,
    col = "blue" # The color of the line when no indiviual letters can be plotted due to size limitations.
  )

  # Set display parameters for the fusion reads
  Gviz::displayPars(fusion@fusionReadsAlignment) <- list(
    showTitle = FALSE,
    min.height = 20, # Minimum height of nucleotide names in pixels
    min.width = 10, # Minimum width of nucleotide names in pixels
    showMismatches = FALSE # Show mismatched reads?
  )

  # Do some plotting

  if (showAllNucleotides) {
    # Create the grid we're gonna use for the plot
    nf <- grid::grid.layout(
      nrow = 1,
      ncol = 1)

    # Open a new page to plot it
    grid::grid.newpage()

    # Divide the page according to the grid we created
    grid::pushViewport(grid::viewport(layout = nf))

    # Open row 1, column 1
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    # Plot the fusion read alignment
    Gviz::plotTracks(
      c(fusion@fusionReadsAlignment, ref_chimera),
      from = 1,
      to = nchar(fusionSequence),
      chromosome = fusion@fusionReadsAlignment@chromosome,
      type = "pileup",
      sizes = c(6, 1),
      add = TRUE
    )
    # Calculate where in the alignment the breakpoint is
    breakpoint <- nchar(fusion@geneA@junctionSequence)/nchar(fusionSequence)
    # Plot line indicating fusion breakpoint
    grid::grid.lines(x = grid::unit(breakpoint, "npc"),
               y = grid::unit(c(0, 1), "npc"),
               default.units = "npc",
               arrow = NULL,
               name = NULL,
               gp = grid::gpar(
                 col = "black",
                 lty = "dotted"),
               draw = TRUE,
               vp = NULL)
    # Close row 1, column 1
    grid::popViewport(1)
  } else {
    # Create the grid we're gonna use for the plot
    nf <- grid::grid.layout(
      nrow = 2,
      ncol = 1,
      heights = grid::unit(c(6, 1), "null"),
      widths = grid::unit(c(1, 1), "null"))

    # Open a new page to plot it
    grid::grid.newpage()

    # Divide the page according to the grid we created
    grid::pushViewport(grid::viewport(layout = nf))

    # Open row 1, column 1
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    # Plot the fusion read alignment
    Gviz::plotTracks(
      fusion@fusionReadsAlignment,
      from = 1,
      to = nchar(fusionSequence),
      chromosome = fusion@fusionReadsAlignment@chromosome,
      type = "pileup",
      add = TRUE
    )
    # Calculate where in the alignment the breakpoint is
    breakpoint <- nchar(fusion@geneA@junctionSequence)/nchar(fusionSequence)
    # Plot line indicating fusion breakpoint
    grid::grid.lines(x = grid::unit(breakpoint, "npc"),
               y = grid::unit(c(0, 1), "npc"),
               default.units = "npc",
               arrow = NULL,
               name = NULL,
               gp = grid::gpar(
                 col = "black",
                 lty = "dotted"),
               draw = TRUE,
               vp = NULL)
    # Close row 1, column 1
    grid::popViewport(1)

    # Open row 1, column 2
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    # Plot a part of the fusion junction sequence
    Gviz::plotTracks(
      ref_chimera_short,
      from = 1,
      to = nchar(fusionSequenceShort),
      chromosome = fusion@fusionReadsAlignment@chromosome,
      add = TRUE
    )
    # Plot line indicating fusion breakpoint
    grid::grid.lines(x = grid::unit(0.5, "npc"),
               y = grid::unit(c(0, 1), "npc"),
               default.units = "npc",
               arrow = NULL,
               name = NULL,
               gp = grid::gpar(
                 col = "black",
                 lty = "dotted"),
               draw = TRUE,
               vp = NULL)
    # Plot line between lines
    grid::grid.lines(x = grid::unit(c(0.5, breakpoint), "npc"),
               y = grid::unit(1, "npc"),
               default.units = "npc",
               arrow = NULL,
               name = NULL,
               gp = grid::gpar(
                 col = "black",
                 lty = "dotted"),
               draw = TRUE,
               vp = NULL)
    # Close row 1, column 2
    grid::popViewport(1)
  }

}
