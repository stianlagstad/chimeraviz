#' Create a plot of the reads supporting the given fusion.
#'
#' This function takes a Fusion object and plots the reads supporting the fusion
#' on top of the fusion sequence (fusion@junction_sequence), provided that
#' add_fusion_reads_alignment() has been run earlier in order to add fusion reads
#' alignment data to the fusion object.
#'
#' Note that the package used for plotting, Gviz, is strict on chromosome names.
#' If the plot produced doesn't show the reads, the problem might be solved by
#' naming the fusion sequence "chrNA".
#'
#' @param fusion The Fusion object to plot.
#' @param show_all_nucleotides Boolean indicating whether or not to show all
#' nucleotides. If FALSE, then only nucleotide_amount amount of nucleotides will
#' be shown on each end of the fusion junction. If TRUE, then the whole fusion
#' junction sequence will be shown.
#' @param nucleotide_amount The number of nucleotides to show on each end of the
#' fusion junction sequence. Defaults to 10. Only applicable if
#' show_all_nucleotides is set to TRUE.
#'
#' @return Creates a fusion reads plot.
#'
#' @seealso add_fusion_reads_alignment
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' # Find the specific fusion we have aligned reads for
#' fusion <- get_fusion_by_id(fusions, 5267)
#' bamfile <- system.file(
#'   "extdata",
#'   "5267readsAligned.bam",
#'   package="chimeraviz")
#' # Add the bam file of aligned fusion reads to the fusion object
#' fusion <- add_fusion_reads_alignment(fusion, bamfile)
#' # Temporary file to store the plot
#' pngFilename <- tempfile(
#'   pattern = "fusionPlot",
#'   fileext = ".png",
#'   tmpdir = tempdir())
#' # Calculate image size based on supporting reads and lenght of junction
#' # sequence.
#' imageWidth <- (nchar(partner_gene_junction_sequence(upstream_partner_gene(fusion))) +
#'   nchar(partner_gene_junction_sequence(downstream_partner_gene(fusion)))) * 15
#' imageHeight <- (fusion_split_reads_count(fusion)+fusion_spanning_reads_count(fusion)) * 20
#' # Open device
#' png(pngFilename, width = imageWidth, height = imageHeight)
#' # Now we can plot
#' plot_fusion_reads(fusion)
#' # Close device
#' dev.off()
#'
#' @export
plot_fusion_reads <- function(
  fusion,
  show_all_nucleotides = TRUE,
  nucleotide_amount = 10) {

  .validate_plot_fusion_reads_params(
    fusion,
    show_all_nucleotides,
    nucleotide_amount
  )

  # Create Biostrings::DNAStringSet of fusion sequence
  fusion_sequence <- Biostrings::DNAStringSet(
    x = c(fusion@gene_upstream@junction_sequence,
          fusion@gene_downstream@junction_sequence))
  names(fusion_sequence) <- "chrNA"
  # And also of a "zoomed-in" version only displaying some nucleotides on each
  # end of the fusion junction
  fusion_sequence_short <- Biostrings::DNAStringSet(
    x = c(
      substr(
        fusion@gene_upstream@junction_sequence,
        nchar(fusion@gene_upstream@junction_sequence) - nucleotide_amount + 1,
        nchar(fusion@gene_upstream@junction_sequence)),
      substr(
        fusion@gene_downstream@junction_sequence,
        nchar(fusion@gene_downstream@junction_sequence) -
          nucleotide_amount + 1 - 1,
        # Added -1 again on the line above because Gviz does not show the last
        # nucleotide. Ref:
        # http://permalink.gmane.org/gmane.science.biology.informatics.conductor.devel/5726 # Exclude Linting
        # Or if that link is down, google "Gviz doesn't provide enough
        # horizontial space"
        nchar(fusion@gene_downstream@junction_sequence))))
  names(fusion_sequence_short) <- "chrNA"

  # Load fusion sequences to tracks
  ref_chimera_short <- Gviz::SequenceTrack(
    fusion_sequence_short,
    genome = fusion@genome_version,
    name = "Fusion Sequence")
  ref_chimera <- Gviz::SequenceTrack(
    fusion_sequence,
    genome = fusion@genome_version,
    name = "Fusion Sequence")

  # DISPLAY PARAMETERS

  # Set display parameters for the fusion sequence
  display_parameters <- list(
    showTitle = FALSE,
    # The color of the line when no indiviual letters can be plotted due to
    # size limitations.
    col = "blue"
  )
  Gviz::displayPars(ref_chimera) <- display_parameters
  Gviz::displayPars(ref_chimera_short) <- display_parameters

  # Set display parameters for the fusion reads
  Gviz::displayPars(fusion@fusion_reads_alignment) <- list(
    showTitle = FALSE,
    showMismatches = FALSE # Show mismatched reads?
  )

  # Do some plotting

  if (show_all_nucleotides) {
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
      c(fusion@fusion_reads_alignment, ref_chimera),
      from = 1,
      to = nchar(fusion_sequence),
      chromosome = fusion@fusion_reads_alignment@chromosome,
      type = "pileup",
      sizes = c(6, 1),
      add = TRUE
    )
    # Calculate where in the alignment the breakpoint is
    breakpoint <-
      nchar(fusion@gene_upstream@junction_sequence) / nchar(fusion_sequence)
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
      fusion@fusion_reads_alignment,
      from = 1,
      to = nchar(fusion_sequence),
      chromosome = fusion@fusion_reads_alignment@chromosome,
      type = "pileup",
      add = TRUE
    )
    # Calculate where in the alignment the breakpoint is
    breakpoint <-
      nchar(fusion@gene_upstream@junction_sequence) / nchar(fusion_sequence)
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

    # Open row 2, column 1
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    # Plot a part of the fusion junction sequence
    Gviz::plotTracks(
      ref_chimera_short,
      from = 1,
      to = nchar(fusion_sequence_short),
      chromosome = fusion@fusion_reads_alignment@chromosome,
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

.validate_plot_fusion_reads_params <- function(
  fusion,
  show_all_nucleotides,
  nucleotide_amount
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- checkmate::makeAssertCollection()

  .is_fusion_valid(argument_checker, fusion)

  # Check if the Gviz::AlignmentsTrack object in
  # fusion@fusion_reads_alignment actually holds data, and is not just the empty
  # object it is after loading the initial fusion data. In other words, make
  # sure add_fusion_reads_alignment() has been run.
  if (class(fusion@fusion_reads_alignment) != "ReferenceAlignmentsTrack") {
    argument_checker$push(
      paste("fusion@fusion_reads_alignment should hold a",
                  "Gviz::ReferenceAlignmentsTrack object.",
                  "See add_fusion_reads_alignment().")
    )
  }

  .is_parameter_boolean(
    argument_checker,
    show_all_nucleotides,
    "show_all_nucleotides")

  .is_nucleotide_amount_valid(
    argument_checker,
    nucleotide_amount,
    fusion
  )

  # Return errors and warnings (if any)
  checkmate::reportAssertions(argument_checker)
}
