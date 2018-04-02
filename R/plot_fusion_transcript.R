#' Plot possible fusion transcripts based on annotation.
#'
#' This function takes a fusion object and an ensembldb object and plots the
#' reduced version of the fusion transcript. This transcript consist of the
#' "mashed together" version of all possible fusion transcripts based on known
#' annotations. If a bamfile is specified, the fusion transcript will be
#' plotted with coverage information.
#'
#' Note that the transcript database used (the edb object) must have the same
#' seqnames as any bamfile used. Otherwise the coverage data will be wrong.
#'
#' @param fusion The Fusion object to plot.
#' @param edb The edb object that will be used to fetch data.
#' @param bamfile The bamfile with RNA-seq data.
#' @param which_transcripts This character vector decides which transcripts are
#' to be plotted. Can be "exonBoundary", "withinExon", "withinIntron",
#' "intergenic", or a character vector with specific transcript ids. Default
#' value is "exonBoundary".
#' @param bedgraphfile A bedGraph file to use instead of the bamfile to plot
#' coverage.
#'
#' @return Creates a fusion transcript plot.
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
#' plot_fusion_transcript(
#'   fusion = fusion,
#'   bamfile = bamfile5267,
#'   edb = edb)
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
#' plot_fusion_transcript(
#'   fusion = fusion,
#'   bamfile = bamfile5267,
#'   edb = edb)
#' # Close device
#' dev.off()
#'
#' @export
plot_fusion_transcript <- function(
  fusion,
  edb = NULL,
  bamfile = NULL,
  which_transcripts = "exonBoundary",
  bedgraphfile = NULL) {

  .validate_plot_fusion_transcript_params(
    fusion,
    edb,
    bamfile,
    which_transcripts,
    bedgraphfile
  )
  fusion <- .get_transcripts_if_not_there(fusion, edb)

  # Select which transcripts to use
  transcripts_upstream <-
    select_transcript(fusion@gene_upstream, which_transcripts)
  transcripts_downstream <-
    select_transcript(fusion@gene_downstream, which_transcripts)

  # Check that there's at least one transcript that has the fusion breakpoint
  # within the transcript.
  .check_that_breakpoints_are_within_transcripts(
    fusion,
    transcripts_upstream,
    transcripts_downstream
  )

  # Reduce transcript into one for each gene
  transcripts_upstream <- reduce(unlist(transcripts_upstream))
  transcripts_downstream <- reduce(unlist(transcripts_downstream))
  # Add transcript metadata column to make Gviz group correctly
  mcols(transcripts_upstream)$transcript <- fusion@gene_upstream@name
  mcols(transcripts_downstream)$transcript <- fusion@gene_downstream@name
  # Add symbol metadata column to make Gviz show ids
  mcols(transcripts_upstream)$symbol <- fusion@gene_upstream@name
  mcols(transcripts_downstream)$symbol <- fusion@gene_downstream@name

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

  if (fusion@gene_upstream@strand == "+") {
    fusion_exons_upstream <-
      transcripts_upstream[
        start(ranges(transcripts_upstream)) <
          fusion@gene_upstream@breakpoint &
        end(ranges(transcripts_upstream)) <=
          fusion@gene_upstream@breakpoint
      ]
  } else {
    fusion_exons_upstream <-
      transcripts_upstream[
        start(ranges(transcripts_upstream)) >=
          fusion@gene_upstream@breakpoint &
        end(ranges(transcripts_upstream)) >
          fusion@gene_upstream@breakpoint
      ]
  }
  if (fusion@gene_downstream@strand == "+") {
    fusion_exons_downstream <-
      transcripts_downstream[
        start(ranges(transcripts_downstream)) >=
          fusion@gene_downstream@breakpoint &
        end(ranges(transcripts_downstream)) >
          fusion@gene_downstream@breakpoint
      ]
  } else {
    fusion_exons_downstream <-
      transcripts_downstream[
        start(ranges(transcripts_downstream)) <
          fusion@gene_downstream@breakpoint &
        end(ranges(transcripts_downstream)) <=
          fusion@gene_downstream@breakpoint
      ]
  }

  # Reverse order of exons if on minus strand
  if (fusion@gene_upstream@strand == "-") {
    fusion_exons_upstream <- rev(fusion_exons_upstream)
  }
  if (fusion@gene_downstream@strand == "-") {
    fusion_exons_downstream <- rev(fusion_exons_downstream)
  }

  # Build fusion transcript
  fusion_transcript <- append(
    fusion_exons_upstream,
    fusion_exons_downstream
  )

  # If we've got a bamfile, calculate coverage
  if (!is.null(bamfile)) {
    # Get coverage
    cov <- GenomicAlignments::coverage(
      bamfile,
      param = Rsamtools::ScanBamParam(which = fusion_transcript))
    # Coverage only for my transcript
    cov <- cov[fusion_transcript]
    # Unlist
    cov <- suppressWarnings(unlist(cov))
    # Turn into numeric vector
    cov <- as.numeric(cov)

    # Create coverage track
    d_track <- DataTrack(
      start = 1:length(cov),
      width = 1,
      chromosome = "chrNA",
      genome = "hg19",
      name = "Coverage",
      type = "h",
      col = "orange",
      fill = "orange",
      data = cov)
    # Set display parameters
    Gviz::displayPars(d_track) <- list(
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
  } else if (!is.null(bedgraphfile)) {
    # We're getting coverage data from a bedGraph file
    d_track <- DataTrack(
      range = bedgraphfile,
      genome = "hg19",
      chromosome = "chrNA",
      name = "Coverage",
      type = "h",
      col = "orange",
      fill = "orange")
    # Set display parameters
    Gviz::displayPars(d_track) <- list(
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
  for (i in 1:length(fusion_transcript)) {
    if (i == 1) {
      fusion_transcript[i] <-
        IRanges::shift(
          fusion_transcript[i],
          shift = 1 - start(fusion_transcript[i])
        )
    } else {
      shift_down_by <-
        - (start(fusion_transcript[i]) - end(fusion_transcript[i - 1])) + 1
      fusion_transcript[i] <-
        IRanges::shift(
          fusion_transcript[i],
          shift = shift_down_by
        )
    }
  }

  # Create new GRanges so that we can set the right seqname
  gr <- GRanges(seqnames = "chrNA", ranges = ranges(fusion_transcript))
  mcols(gr)$symbol <- mcols(fusion_transcript)$symbol

  # Create transcript track
  tr_track <- Gviz::GeneRegionTrack(
    gr,
    chromosome = "chrNA")
  # Set display parameters
  Gviz::displayPars(tr_track) <- list(
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
  Gviz::feature(tr_track)[
    Gviz::symbol(tr_track) == fusion@gene_upstream@name
  ] <- "includedA"
  Gviz::feature(tr_track)[
    Gviz::symbol(tr_track) == fusion@gene_downstream@name
  ] <- "includedB"
  Gviz::displayPars(tr_track) <- interestcolor

  # Create axis track
  axis_track <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axis_track) <- list(
    col = "black", # line color
    fontcolor = "black",
    fontsize = 8,
    lwd = 1.5
  )

  # If we've got a bamfile, plot with coverage
  if (!is.null(bamfile) || !is.null(bedgraphfile)) {
    Gviz::plotTracks(
      list(d_track, tr_track, axis_track),
      main = paste(
        fusion@gene_upstream@name,
        fusion@gene_downstream@name,
        sep = ":"
      ),
      chromosome = "chrNA")
  } else {
    Gviz::plotTracks(
      list(tr_track, axis_track),
      main = paste(
        fusion@gene_upstream@name,
        fusion@gene_downstream@name,
        sep = ":"
      ),
      chromosome = "chrNA")
  }

}

.validate_plot_fusion_transcript_params <- function(
  fusion,
  edb,
  bamfile,
  which_transcripts,
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

  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)
}
