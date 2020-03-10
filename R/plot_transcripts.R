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
#' @param which_transcripts This character vector decides which transcripts are
#' to be plotted. Can be "exonBoundary", "withinExon", "withinIntron",
#' "intergenic", or a character vector with specific transcript ids. Default
#' value is "exonBoundary".
#' @param ylim Limits for the coverage y-axis.
#' @param non_ucsc Boolean indicating whether or not the bamfile used has UCSC-
#' styled chromosome names (i.e. with the "chr" prefix). Setting this to true
#' lets you use a bamfile with chromosome names like "1" and "X", instead of
#' "chr1" and "chrX".
#' @param reduce_transcripts Boolean indicating whether or not to reduce all
#' transcripts into a single transcript for each partner gene.
#' @param bedgraphfile A bedGraph file to use instead of the bamfile to plot
#' coverage.
#'
#' @return Creates a fusion transcripts plot.
#'
#' @examples
#' # Load data and example fusion event
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
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
#' plot_transcripts(
#'   fusion = fusion,
#'   edb = edb,
#'   bamfile = bamfile5267,
#'   non_ucsc = TRUE)
#' # Close device
#' dev.off()
#'
#' # Example using a .bedGraph file instead of a .bam file:
#' # Load data and example fusion event
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Load edb
#' edbSqliteFile <- system.file(
#'   "extdata",
#'   "Homo_sapiens.GRCh37.74.sqlite",
#'   package="chimeraviz")
#' edb <- ensembldb::EnsDb(edbSqliteFile)
#' # bedgraphfile with coverage data from the regions of this fusion event
#' bedgraphfile <- system.file(
#'   "extdata",
#'   "fusion5267and11759reads.bedGraph",
#'   package="chimeraviz")
#' # Temporary file to store the plot
#' pngFilename <- tempfile(
#'   pattern = "fusionPlot",
#'   fileext = ".png",
#'   tmpdir = tempdir())
#' # Open device
#' png(pngFilename, width = 500, height = 500)
#' # Plot!
#' plot_transcripts(
#'   fusion = fusion,
#'   edb = edb,
#'   bedgraphfile = bedgraphfile,
#'   non_ucsc = TRUE)
#' # Close device
#' dev.off()
#'
#' @importFrom ensembldb EnsDb
#'
#' @export
plot_transcripts <- function(
  fusion,
  edb = NULL,
  bamfile = NULL,
  which_transcripts = "exonBoundary",
  non_ucsc = TRUE,
  ylim = c(0, 1000),
  reduce_transcripts = FALSE,
  bedgraphfile = NULL) {

  .validate_plot_fusion_transcripts_params(
    fusion,
    edb,
    bamfile,
    which_transcripts,
    non_ucsc,
    ylim,
    reduce_transcripts,
    bedgraphfile)

  fusion <- .get_transcripts_if_not_there(fusion, edb)

  # Select which transcripts to use
  transcripts_upstream <-
    select_transcript(fusion@gene_upstream, which_transcripts)
  transcripts_downstream <-
    select_transcript(fusion@gene_downstream, which_transcripts)

  # Set seqlevels to UCSC
  GenomeInfoDb::seqlevelsStyle(transcripts_upstream) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(transcripts_downstream) <- "UCSC"

  if (reduce_transcripts) {
    reduced_transcripts_upstream <-
      GenomicRanges::reduce(unlist(transcripts_upstream))
    reduced_transcripts_downstream <-
      GenomicRanges::reduce(unlist(transcripts_downstream))
    # Add transcript metadata column to make Gviz group correctly
    mcols(reduced_transcripts_upstream)$transcript <-
      fusion@gene_upstream@name
    mcols(reduced_transcripts_downstream)$transcript <-
      fusion@gene_downstream@name
    # Create the tracks
    tr_upstream_track <- Gviz::GeneRegionTrack(
      reduced_transcripts_upstream,
      chromosome = fusion@gene_upstream@chromosome
    )
    tr_downstream_track <- Gviz::GeneRegionTrack(
      reduced_transcripts_downstream,
      chromosome = fusion@gene_downstream@chromosome
    )
    # Set symbol names on the tracks
    symbol(tr_upstream_track) <- fusion@gene_upstream@name
    symbol(tr_downstream_track) <- fusion@gene_downstream@name
  } else {
    # Create the tracks
    tr_upstream_track <- Gviz::GeneRegionTrack(
      data.frame(transcripts_upstream),
      chromosome = fusion@gene_upstream@chromosome)
    tr_downstream_track <- Gviz::GeneRegionTrack(
      data.frame(transcripts_downstream),
      chromosome = fusion@gene_downstream@chromosome)
  }

  # Set display parameters
  transcript_display_params <- list(
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
  Gviz::displayPars(tr_upstream_track) <- transcript_display_params
  Gviz::displayPars(tr_downstream_track) <- transcript_display_params

  # Create axis track
  axis_track <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axis_track) <- list(
    labelPos = "below",
    cex = 1,
    col = "black", # line color
    fontcolor = "black",
    lwd = 1.5
  )

  # Helper variables
  gene_upstream_at_minus_strand <- fusion@gene_upstream@strand == "-"
  gene_downstream_at_minus_strand <- fusion@gene_downstream@strand == "-"

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
  if (reduce_transcripts) {
    Gviz::feature(tr_upstream_track) <- "protein_coding_excludedA"
    Gviz::feature(tr_downstream_track) <- "protein_coding_excludedB"
  } else {
    Gviz::feature(tr_upstream_track)[
      Gviz::feature(tr_upstream_track) == "5utr"
    ] <- "5utr_excludedA"
    Gviz::feature(tr_upstream_track)[
      Gviz::feature(tr_upstream_track) == "protein_coding"
    ] <- "protein_coding_excludedA"
    Gviz::feature(tr_upstream_track)[
      Gviz::feature(tr_upstream_track) == "3utr"
    ] <- "3utr_excludedA"

    Gviz::feature(tr_downstream_track)[
      Gviz::feature(tr_downstream_track) == "5utr"
    ] <- "5utr_excludedB"
    Gviz::feature(tr_downstream_track)[
      Gviz::feature(tr_downstream_track) == "protein_coding"
    ] <- "protein_coding_excludedB"
    Gviz::feature(tr_downstream_track)[
      Gviz::feature(tr_downstream_track) == "3utr"
    ] <- "3utr_excludedB"
  }

  if (gene_upstream_at_minus_strand) {
    message(paste("As ",
                  fusion@gene_upstream@name,
                  " is on the minus strand, the plot for this gene will be ",
                  "reversed",
                  sep = ""))
    highlight_start_upstream <- fusion@gene_upstream@breakpoint
    highlight_end_upstream <-
      max(start(tr_upstream_track) + width(tr_upstream_track) - 1)

    # Color the exons included/not included in the fusion
    Gviz::feature(tr_upstream_track)[
      start(tr_upstream_track) <= highlight_end_upstream &
        end(tr_upstream_track) >= highlight_start_upstream &
        Gviz::feature(tr_upstream_track) == "5utr_excludedA"
    ] <- "5utr_includedA"

    Gviz::feature(tr_upstream_track)[
      start(tr_upstream_track) <= highlight_end_upstream &
        end(tr_upstream_track) >= highlight_start_upstream &
        Gviz::feature(tr_upstream_track) == "protein_coding_excludedA"
    ] <- "protein_coding_includedA"

    Gviz::feature(tr_upstream_track)[
      start(tr_upstream_track) <= highlight_end_upstream &
        end(tr_upstream_track) >= highlight_start_upstream &
        Gviz::feature(tr_upstream_track) == "3utr_excludedA"
    ] <- "3utr_includedA"

  } else {
    highlight_start_upstream <- min(start(tr_upstream_track))
    highlight_end_upstream <- fusion@gene_upstream@breakpoint

    # Color the exons included/not included in the fusion
    Gviz::feature(tr_upstream_track)[
      start(tr_upstream_track) >= highlight_start_upstream &
      end(tr_upstream_track) <= highlight_end_upstream &
      Gviz::feature(tr_upstream_track) == "5utr_excludedA"
    ] <- "5utr_includedA"

    Gviz::feature(tr_upstream_track)[
      start(tr_upstream_track) >= highlight_start_upstream &
      end(tr_upstream_track) <= highlight_end_upstream &
      Gviz::feature(tr_upstream_track) == "protein_coding_excludedA"
    ] <- "protein_coding_includedA"

    Gviz::feature(tr_upstream_track)[
      start(tr_upstream_track) >= highlight_start_upstream &
      end(tr_upstream_track) <= highlight_end_upstream &
      Gviz::feature(tr_upstream_track) == "3utr_excludedA"
    ] <- "3utr_includedA"

  }

  if (gene_downstream_at_minus_strand) {
    message(paste("As ",
                  fusion@gene_downstream@name,
                  " is on the minus strand, the plot for this gene will be ",
                  "reversed",
                  sep = ""))
    highlight_start_downstream <- min(start(tr_downstream_track))
    highlight_end_downstream <- fusion@gene_downstream@breakpoint

    # Color the exons included/not included in the fusion
    Gviz::feature(tr_downstream_track)[
      start(tr_downstream_track) <= highlight_end_downstream &
      end(tr_downstream_track) >= highlight_start_downstream &
      Gviz::feature(tr_downstream_track) == "5utr_excludedB"
    ] <- "5utr_includedB"

    Gviz::feature(tr_downstream_track)[
      start(tr_downstream_track) <= highlight_end_downstream &
      end(tr_downstream_track) >= highlight_start_downstream &
      Gviz::feature(tr_downstream_track) == "protein_coding_excludedB"
    ] <- "protein_coding_includedB"

    Gviz::feature(tr_downstream_track)[
      start(tr_downstream_track) <= highlight_end_downstream &
      end(tr_downstream_track) >= highlight_start_downstream &
      Gviz::feature(tr_downstream_track) == "3utr_excludedB"
    ] <- "3utr_includedB"

  } else {
    highlight_start_downstream <- fusion@gene_downstream@breakpoint
    highlight_end_downstream <-
      max(start(tr_downstream_track) + width(tr_downstream_track) - 1)

    # Color the exons included/not included in the fusion
    Gviz::feature(tr_downstream_track)[
      start(tr_downstream_track) >= highlight_start_downstream &
      end(tr_downstream_track) <= highlight_end_downstream &
      Gviz::feature(tr_downstream_track) == "5utr_excludedB"
    ] <- "5utr_includedB"

    Gviz::feature(tr_downstream_track)[
      start(tr_downstream_track) >= highlight_start_downstream &
      end(tr_downstream_track) <= highlight_end_downstream &
      Gviz::feature(tr_downstream_track) == "protein_coding_excludedB"
    ] <- "protein_coding_includedB"

    Gviz::feature(tr_downstream_track)[
      start(tr_downstream_track) >= highlight_start_downstream &
      end(tr_downstream_track) <= highlight_end_downstream &
      Gviz::feature(tr_downstream_track) == "3utr_excludedB"
    ] <- "3utr_includedB"

  }

  # Update transcript track with colors
  Gviz::displayPars(tr_upstream_track) <- interestcolor
  Gviz::displayPars(tr_downstream_track) <- interestcolor

  # Highlight transcript tracks

  # Display parameters for highlight tracks
  display_pars_highlight_tracks <- list(
    fill = "lightgrey",
    col = "lightgrey"
  )

  # Only plot highlight tracks if we actually have transcripts that has the
  # fusion breakpoint within them

  if (any(start(tr_upstream_track) < fusion@gene_upstream@breakpoint) &&
      any(fusion@gene_upstream@breakpoint < end(tr_upstream_track))) {

    gr_track_highlight_upstream <- Gviz::HighlightTrack(
      trackList = tr_upstream_track,
      start = highlight_start_upstream,
      end = highlight_end_upstream,
      chromosome = fusion@gene_upstream@chromosome)
    Gviz::displayPars(gr_track_highlight_upstream) <-
      display_pars_highlight_tracks
  }

  if (any(start(tr_downstream_track) < fusion@gene_downstream@breakpoint) &&
      any(fusion@gene_downstream@breakpoint < end(tr_downstream_track))) {

    gr_track_highlight_downstream <- Gviz::HighlightTrack(
      trackList = tr_downstream_track,
      start = highlight_start_downstream,
      end = highlight_end_downstream,
      chromosome = fusion@gene_downstream@chromosome)
    Gviz::displayPars(gr_track_highlight_downstream) <-
      display_pars_highlight_tracks
  }

  # Plot coverage?
  if (!is.null(bamfile) || !is.null(bedgraphfile)) {
    # Create alignment track
    if (!is.null(bamfile)) {
      # We're getting coverage data from a bam file
      if (non_ucsc) {
        # If the bam file has non-ucsc chromosome names, i.e. "1" instead of
        # "chr1", then we need to use a custom import function that adds "chr"
        # to the chromosome names, to keep Gviz happy. Gviz strongly prefers
        # having the "chr" prefix.
        al_track <- Gviz::AlignmentsTrack(
          bamfile,
          isPaired = TRUE,
          genome = fusion@genome_version,
          ylim = ylim,
          importFunction = import_function_non_ucsc)
      } else {
        al_track <- Gviz::AlignmentsTrack(
          bamfile,
          isPaired = TRUE,
          genome = fusion@genome_version,
          ylim = ylim)
      }
      # Set display paramters
      Gviz::displayPars(al_track) <- list(
        showTitle = FALSE, # hide name of track
        background.panel = "transparent", # bg color of the content panel
        background.title = "transparent", # bg color for the title panels
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
    } else {
      # We're getting coverage data from a bedGraph file
      al_track <- DataTrack(
        range = bedgraphfile,
        ylim = ylim,
        genome = "hg19",
        name = "Coverage",
        type = "h",
        col = "orange",
        fill = "orange")
      # Set display parameters
      Gviz::displayPars(al_track) <- list(
        showTitle = FALSE, # hide name of track
        background.panel = "transparent", # bg color of the content panel
        background.title = "transparent", # bg color for the title panels
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

    # Highlight coverage tracks

    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(tr_upstream_track) < fusion@gene_upstream@breakpoint) &&
        any(fusion@gene_upstream@breakpoint < end(tr_upstream_track))) {

      al_track_highlight_upstream <- Gviz::HighlightTrack(
        trackList = al_track,
        start = highlight_start_upstream,
        end = highlight_end_upstream,
        chromosome = fusion@gene_upstream@chromosome)
      Gviz::displayPars(al_track_highlight_upstream) <-
        display_pars_highlight_tracks
    }

    if (any(start(tr_downstream_track) < fusion@gene_downstream@breakpoint) &&
        any(fusion@gene_downstream@breakpoint < end(tr_downstream_track))) {

      al_track_highlight_downstream <- Gviz::HighlightTrack(
        trackList = al_track,
        start = highlight_start_downstream,
        end = highlight_end_downstream,
        chromosome = fusion@gene_downstream@chromosome)
      Gviz::displayPars(al_track_highlight_downstream) <-
        display_pars_highlight_tracks
    }
  }

  # Do some plotting

  # Calculate heights in the plot grid according to how many transcripts there
  # are in each partnergene
  if (reduce_transcripts) {
    number_of_transcripts_upstream <- 1
    number_of_transcripts_downstream <- 1
  } else {
    number_of_transcripts_upstream <- length(transcripts_upstream)
    number_of_transcripts_downstream <- length(transcripts_downstream)
  }

  # Create the grid we're gonna use for the plot
  nf <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(
      c(
        number_of_transcripts_upstream,
        number_of_transcripts_downstream
      ),
      "null"
    )
  )

  # Open a new page to plot it
  grid::grid.newpage()

  # Divide the page according to the grid we created
  grid::pushViewport(grid::viewport(layout = nf))

  # Open row 1, column 1
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  # Plot transcriptsA and coverage with highlight
  # Plot coverage?
  if (!is.null(bamfile) || !is.null(bedgraphfile)) {

    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(tr_upstream_track) < fusion@gene_upstream@breakpoint) &&
        any(fusion@gene_upstream@breakpoint < end(tr_upstream_track))) {

      Gviz::plotTracks(
        # without this gviz create cluster_X entries in the GeneRegionTrack
        collapse = FALSE,
        list(gr_track_highlight_upstream, al_track_highlight_upstream),
        from = min(start(tr_upstream_track)),
        to = max(end(tr_upstream_track)),
        sizes = c(5, 2),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = gene_upstream_at_minus_strand)
    } else {
      Gviz::plotTracks(
        # without this gviz create cluster_X entries in the GeneRegionTrack
        collapse = FALSE,
        list(tr_upstream_track, al_track),
        from = min(start(tr_upstream_track)),
        to = max(end(tr_upstream_track)),
        sizes = c(5, 2),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = gene_upstream_at_minus_strand)
    }
  } else {
    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(tr_upstream_track) < fusion@gene_upstream@breakpoint) &&
        any(fusion@gene_upstream@breakpoint < end(tr_upstream_track))) {

      Gviz::plotTracks(
        # without this gviz create cluster_X entries in the GeneRegionTrack
        collapse = FALSE,
        gr_track_highlight_upstream,
        from = min(start(tr_upstream_track)),
        to = max(end(tr_upstream_track)),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = gene_upstream_at_minus_strand)
    } else {
      Gviz::plotTracks(
        # without this gviz create cluster_X entries in the GeneRegionTrack
        collapse = FALSE,
        tr_upstream_track,
        from = min(start(tr_upstream_track)),
        to = max(end(tr_upstream_track)),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = gene_upstream_at_minus_strand)
    }
  }
  # Close row 1, column 1
  grid::popViewport(1)

  # Open row 2, column 1
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  # Plot transcriptsB and coverage with highlight
  # Plot coverage?
  if (!is.null(bamfile) || !is.null(bedgraphfile)) {
    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(tr_downstream_track) < fusion@gene_downstream@breakpoint) &&
        any(fusion@gene_downstream@breakpoint < end(tr_downstream_track))) {

      Gviz::plotTracks(
        # without this gviz create cluster_X entries in the GeneRegionTrack
        collapse = FALSE,
        list(gr_track_highlight_downstream, al_track_highlight_downstream),
        from = min(start(tr_downstream_track)),
        to = max(end(tr_downstream_track)),
        sizes = c(5, 2),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = gene_downstream_at_minus_strand)
    } else {
      Gviz::plotTracks(
        # without this gviz create cluster_X entries in the GeneRegionTrack
        collapse = FALSE,
        list(tr_downstream_track, al_track),
        from = min(start(tr_downstream_track)),
        to = max(end(tr_downstream_track)),
        sizes = c(5, 2),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = gene_downstream_at_minus_strand)
    }
  } else {
    # Only plot highlight tracks if we actually have transcripts that has the
    # fusion breakpoint within them
    if (any(start(tr_downstream_track) < fusion@gene_downstream@breakpoint) &&
        any(fusion@gene_downstream@breakpoint < end(tr_downstream_track))) {

      Gviz::plotTracks(
        # without this gviz create cluster_X entries in the GeneRegionTrack
        collapse = FALSE,
        gr_track_highlight_downstream,
        from = min(start(tr_downstream_track)),
        to = max(end(tr_downstream_track)),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = gene_downstream_at_minus_strand)
    } else {
      Gviz::plotTracks(
        # without this gviz create cluster_X entries in the GeneRegionTrack
        collapse = FALSE,
        tr_downstream_track,
        from = min(start(tr_downstream_track)),
        to = max(end(tr_downstream_track)),
        add = TRUE,
        margin = 3,
        innerMargin = 0,
        # Plot reverse if gene is at minus strand
        reverseStrand = gene_downstream_at_minus_strand)
    }
  }
  # Close row 1, column 1
  grid::popViewport(1)
}

.validate_plot_fusion_transcripts_params <- function(
  fusion,
  edb,
  bamfile,
  which_transcripts,
  non_ucsc,
  ylim,
  reduce_transcripts,
  bedgraphfile
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()

  # Check parameters
  argument_checker <- .is_fusion_valid(argument_checker, fusion)
  argument_checker <- .is_edb_valid(argument_checker, edb, fusion)
  # Either bamfile or bedgraphfile can be given, not both
  bamfile_given <- !is.null(bamfile)
  bedgraphfile_given <- !is.null(bedgraphfile)
  if (bamfile_given && bedgraphfile_given) {
    ArgumentCheck::addError(
      msg = "Either 'bamfile' or 'bedgraphfile' must be given, not both.",
      argcheck = argument_checker
    )
  }
  if (bamfile_given) {
    argument_checker <- .is_bamfile_valid(argument_checker, bamfile)
  }
  if (bedgraphfile_given) {
    argument_checker <- .is_bedgraphfile_valid(argument_checker, bedgraphfile)
  }
  argument_checker <- .is_which_transcripts_valid(
    argument_checker,
    which_transcripts,
    fusion)
  argument_checker <- .is_ylim_valid(argument_checker, ylim)
  argument_checker <- .is_parameter_boolean(
    argument_checker,
    non_ucsc,
    "non_ucsc")
  argument_checker <- .is_parameter_boolean(
    argument_checker,
    reduce_transcripts,
    "reduce_transcripts")

  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)
}
