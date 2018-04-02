#' Plot a fusion event with transcripts, coverage and ideograms.
#'
#' This function creates a plot with information about transcripts, coverage,
#' location and more.
#'
#' plot_fusion() will dispatch to either plot_fusion_separate() or
#' plot_fusion_together(). plot_fusion_separate() will plot the fusion gene
#' partners in separate graphs shown next to each other, while
#' plot_fusion_together() will plot the fusion gene partners in the same graph
#' with the same x-axis. plot_fusion() will dispatch to plot_fusion_together() if
#' the fusion gene partners are on the same strand, same chromosome and are
#' close together (<=50,000 bp apart).
#'
#' @param fusion The Fusion object to plot.
#' @param edb The ensembldb object that will be used to fetch data.
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
#' @return Creates a fusion plot.
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
#' png(pngFilename, width = 1000, height = 750)
#' # Plot!
#' plot_fusion(
#'   fusion = fusion,
#'   bamfile = bamfile5267,
#'   edb = edb,
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
#' png(pngFilename, width = 1000, height = 750)
#' # Plot!
#' plot_fusion(
#'   fusion = fusion,
#'   bedgraphfile = bedgraphfile,
#'   edb = edb,
#'   non_ucsc = TRUE)
#' # Close device
#' dev.off()
#'
#' @export
plot_fusion <- function(
  fusion,
  edb = NULL,
  bamfile = NULL,
  which_transcripts = "exonBoundary",
  ylim = c(0, 1000),
  non_ucsc = TRUE,
  reduce_transcripts = FALSE,
  bedgraphfile = NULL) {

  .validate_plot_fusion_params(
    fusion,
    edb,
    bamfile,
    which_transcripts,
    ylim,
    non_ucsc,
    reduce_transcripts,
    bedgraphfile
  )
  fusion <- .get_transcripts_if_not_there(fusion, edb)

  # Decide whether or not to plot the fusion on the same x-axis
  if (fusion@gene_upstream@chromosome == fusion@gene_downstream@chromosome) {
    message("Fusion is intrachromosomal..")
    if (
      as.character(fusion@gene_upstream@strand) !=
      as.character(fusion@gene_downstream@strand)
    ) {
      message("..but on different strands. Plot separate!")
      plot_fusion_separate(
        fusion = fusion,
        edb = edb,
        bamfile = bamfile,
        which_transcripts = which_transcripts,
        ylim = ylim,
        non_ucsc = non_ucsc,
        reduce_transcripts = reduce_transcripts,
        bedgraphfile = bedgraphfile)
    } else if (
      abs(
        fusion@gene_upstream@breakpoint - fusion@gene_downstream@breakpoint
      ) > 50000
    ) {
      message("..but too far apart. Plot separate!")
      plot_fusion_separate(
        fusion = fusion,
        edb = edb,
        bamfile = bamfile,
        which_transcripts = which_transcripts,
        ylim = ylim,
        non_ucsc = non_ucsc,
        reduce_transcripts = reduce_transcripts,
        bedgraphfile = bedgraphfile)
    } else {
      message("..and close together..")
      # If both genes are on the minus strand, and the location of geneA is to
      # the right of geneB, then we still want to plot them separatly, because
      # we want geneA (the downstream fusion gene partner) to be on the left
      # side of the plot.
      if (fusion@gene_upstream@strand == "-" &&
          fusion@gene_downstream@strand == "-" &&
          min(start(unlist(fusion@gene_upstream@transcripts))) <
          min(start(unlist(fusion@gene_downstream@transcripts)))) {
        message(
          paste0(
            "..but since both genes are on the minus strand the x-axis will ",
            "be flipped in the plot. And since geneA will then be located to ",
            "\" the right of\" geneB, we will still plot this separatly. Plot",
            " separate!"
          )
        )
        plot_fusion_separate(
          fusion = fusion,
          edb = edb,
          bamfile = bamfile,
          which_transcripts = which_transcripts,
          ylim = ylim,
          non_ucsc = non_ucsc,
          reduce_transcripts = reduce_transcripts,
          bedgraphfile = bedgraphfile)
      } else {
        message("Plot together!")
        plot_fusion_together(
          fusion = fusion,
          edb = edb,
          bamfile = bamfile,
          which_transcripts = which_transcripts,
          ylim = ylim,
          non_ucsc = non_ucsc,
          reduce_transcripts = reduce_transcripts,
          bedgraphfile = bedgraphfile)
      }
    }
  } else {
    message("Fusion is interchromosomal. Plot separate!")
    plot_fusion_separate(
      fusion = fusion,
      edb = edb,
      bamfile = bamfile,
      which_transcripts = which_transcripts,
      ylim = ylim,
      non_ucsc = non_ucsc,
      reduce_transcripts = reduce_transcripts,
      bedgraphfile = bedgraphfile)
  }
}

#' @name plot_fusion
#'
#' @export
plot_fusion_separate <- function(
  fusion,
  edb,
  bamfile = NULL,
  which_transcripts = "exonBoundary",
  ylim = c(0, 1000),
  non_ucsc = TRUE,
  reduce_transcripts = FALSE,
  bedgraphfile = NULL) {

  .validate_plot_fusion_params(
    fusion,
    edb,
    bamfile,
    which_transcripts,
    ylim,
    non_ucsc,
    reduce_transcripts,
    bedgraphfile
  )
  fusion <- .get_transcripts_if_not_there(fusion, edb)

  # Create tracks

  # Select transcripts
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

  # Set seqlevels to UCSC
  GenomeInfoDb::seqlevelsStyle(transcripts_upstream) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(transcripts_downstream) <- "UCSC"
  if (reduce_transcripts) {
    reduced_transcripts_upstream <- GenomicRanges::reduce(
      unlist(transcripts_upstream)
    )
    reduced_transcripts_downstream <- GenomicRanges::reduce(
      unlist(transcripts_downstream)
    )
    # Add transcript metadata column to make Gviz group correctly
    mcols(reduced_transcripts_upstream)$transcript <-
      fusion@gene_upstream@name
    mcols(reduced_transcripts_downstream)$transcript <-
      fusion@gene_downstream@name
    # Create the tracks
    tr_upstream_track <- Gviz::GeneRegionTrack(
      data.frame(reduced_transcripts_upstream),
      chromosome = fusion@gene_upstream@chromosome)
    tr_downstream_track <- Gviz::GeneRegionTrack(
      data.frame(reduced_transcripts_downstream),
      chromosome = fusion@gene_downstream@chromosome)
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
  Gviz::displayPars(tr_upstream_track) <- transcript_display_params
  Gviz::displayPars(tr_downstream_track) <- transcript_display_params

  # Create axis track
  axis_track <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axis_track) <- list(
    add53 = TRUE,
    add35 = TRUE,
    col = "black", # line color
    fontcolor = "black",
    fontsize = 15
  )

  # Create ideogram tracks

  # Read cytoband information depending on genome version
  if (fusion@genome_version == "hg19") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.HG19.Human.CytoBandIdeogram.txt",
      package = "chimeraviz")
  } else if (fusion@genome_version == "hg38") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.HG38.Human.CytoBandIdeogram.txt",
      package = "chimeraviz")
  } else if (fusion@genome_version == "mm10") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.MM10.Mus.musculus.CytoBandIdeogram.txt",
      package = "chimeraviz")
  } else {
    stop("Unsupported genome version")
  }
  cytoband <- utils::read.table(cytoband_file)
  # Set names to what Gviz expects
  names(cytoband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

  ideogram_track_downstream <- Gviz::IdeogramTrack(
    genome = fusion@genome_version,
    chromosome = fusion@gene_upstream@chromosome,
    bands = cytoband)
  ideogram_track_downstream <- Gviz::IdeogramTrack(
    genome = fusion@genome_version,
    chromosome = fusion@gene_downstream@chromosome,
    bands = cytoband)
  # Set display paramters
  ideogram_display_params <- list(
    showBandId = TRUE,
    showId = TRUE,
    fontcolor = "black",
    cex = 1.5)
  Gviz::displayPars(ideogram_track_downstream) <- ideogram_display_params
  Gviz::displayPars(ideogram_track_downstream) <- ideogram_display_params

  # Create alignment track
  if (!is.null(bamfile)) {
    # We're getting coverage data from a bam file
    if (non_ucsc) {
      # If the bam file has non-ucsc chromosome names, i.e. "1" instead of
      # "chr1", then we need to use a custom import function that adds "chr" to
      # the chromosome names, to keep Gviz happy. Gviz strongly prefers having
      # the "chr" prefix.
      alignment_track <- Gviz::AlignmentsTrack(
        bamfile,
        isPaired = TRUE,
        genome = fusion@genome_version,
        name = "RNA coverage",
        ylim = ylim,
        importFunction = import_function_non_ucsc)
    } else {
      alignment_track <- Gviz::AlignmentsTrack(
        bamfile,
        isPaired = TRUE,
        genome = fusion@genome_version,
        name = "RNA coverage",
        ylim = ylim)
    }
    # Set display paramters
    Gviz::displayPars(alignment_track) <- list(
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
  } else {
    # We're getting coverage data from a bedGraph file
    alignment_track <- DataTrack(
      range = bedgraphfile,
      ylim = ylim,
      genome = "hg19",
      name = "Coverage",
      type = "h",
      col = "orange",
      fill = "orange")
    # Set display parameters
    Gviz::displayPars(alignment_track) <- list(
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

  # Display parameters for highlight tracks
  display_pars_highlight_tracks <- list(
    fill = "lightgrey",
    col = "lightgrey",
    showTitle = FALSE
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

  gr_track_highlight_upstream <- Gviz::HighlightTrack(
    trackList = tr_upstream_track,
    start = highlight_start_upstream,
    end = highlight_end_upstream,
    chromosome = fusion@gene_upstream@chromosome)
  Gviz::displayPars(gr_track_highlight_upstream) <-
    display_pars_highlight_tracks

  gr_track_highlight_downstream <- Gviz::HighlightTrack(
    trackList = tr_downstream_track,
    start = highlight_start_downstream,
    end = highlight_end_downstream,
    chromosome = fusion@gene_downstream@chromosome)
  Gviz::displayPars(gr_track_highlight_downstream) <-
    display_pars_highlight_tracks

  # Highlight coverage tracks

  al_track_highlight_upstream <- Gviz::HighlightTrack(
    trackList = alignment_track,
    start = highlight_start_upstream,
    end = highlight_end_upstream,
    chromosome = fusion@gene_upstream@chromosome)
  Gviz::displayPars(al_track_highlight_upstream) <-
    display_pars_highlight_tracks

  al_track_highlight_downstream <- Gviz::HighlightTrack(
    trackList = alignment_track,
    start = highlight_start_downstream,
    end = highlight_end_downstream,
    chromosome = fusion@gene_downstream@chromosome)
  Gviz::displayPars(al_track_highlight_downstream) <-
    display_pars_highlight_tracks

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
  # we want to show. If reduce_transcripts=TRUE, then we only want to display
  # one transcript per gene.
  if (reduce_transcripts) {
    transcripts_row_height <- 1
  } else {
    # It depends on max(amountOfTranscriptsA, amountOfTranscriptsB)
    transcripts_row_height <- max(
      length(transcripts_upstream),
      length(transcripts_downstream))
    # Set a max value
    if (transcripts_row_height > 10) {
      transcripts_row_height <- 10
    }
  }

  # If the range for the plot is too short, then GenomeAxisTrack will show
  # ticks too close to each other. This doesn't look nice. A way to solve it is
  # to extend the plot range if the transcript is too short:
  plot_range <- max(end(tr_upstream_track)) - min(start(tr_upstream_track))
  if (plot_range < 20000) {
    offset <- 2500
  } else {
    offset <- 0
  }

  # Open row 1, column 1
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  # Plot transcriptsA and coverage with highlight
  res1 <- Gviz::plotTracks(
    # without this gviz create cluster_X entries in the GeneRegionTrack
    collapse = FALSE,
    c(
      ideogram_track_downstream,
      gr_track_highlight_upstream,
      al_track_highlight_upstream,
      axis_track
    ),
    # 10k added so that we can see all transcript names
    from = min(start(tr_upstream_track)) - offset - 10000,
    # 10k added so that we can see all transcript names
    to = max(end(tr_upstream_track)) + offset + 10000,
    sizes = c(1, transcripts_row_height, 2, .5),
    add = TRUE,
    background.title = "transparent",
    chromosome = rep(fusion@gene_upstream@chromosome, 4),
    # Plot reverse if gene is at minus strand
    reverseStrand = gene_upstream_at_minus_strand)
  # Close row 1, column 1
  grid::popViewport(1)

  # If the range for the plot is too short, then GenomeAxisTrack will show
  # ticks too close to each other. This doesn't look nice. A way to solve it is
  # to extend the plot range if the transcript is too short:
  plot_range <- max(end(tr_downstream_track)) - min(start(tr_downstream_track))
  if (plot_range < 20000) {
    offset <- 2500
  } else {
    offset <- 0
  }

  # Open row 1, column 2
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  # Plot transcriptsB and coverage with highlight
  res2 <- Gviz::plotTracks(
    # without this gviz create cluster_X entries in the GeneRegionTrack
    collapse = FALSE,
    c(
      ideogram_track_downstream,
      gr_track_highlight_downstream,
      al_track_highlight_downstream,
      axis_track
    ),
    # 10k added so that we can see all transcript names
    from = min(start(tr_downstream_track)) - offset - 10000,
    # 10k added so that we can see all transcript names
    to = max(end(tr_downstream_track)) + offset + 10000,
    sizes = c(1, transcripts_row_height, 2, .5),
    add = TRUE,
    background.title = "transparent",
    chromosome = rep(fusion@gene_downstream@chromosome, 4),
    # Plot reverse if gene is at minus strand
    reverseStrand = gene_downstream_at_minus_strand)
  # Close row 1, column 2
  grid::popViewport(1)

  # Draw bezier curve if we're plotting exon boundary transcripts

  upstream_boundary_exons <- fusion@gene_upstream@transcripts[
    mcols(fusion@gene_upstream@transcripts)$transcript_category ==
      "exonBoundary"
  ]
  downstream_boundary_exons <- fusion@gene_downstream@transcripts[
    mcols(fusion@gene_downstream@transcripts)$transcript_category ==
      "exonBoundary"
  ]
  if (length(upstream_boundary_exons) > 0 &
      length(downstream_boundary_exons) > 0) {

    # We need to figure out which exons in tr_upstream_track and
    # tr_downstream_track matches the fusion breakpoints
    # fusion@gene_upstream@breakpoint and fusion@gene_downstream@breakpoint:
    if (gene_upstream_at_minus_strand) {
      exon_breakpoint_upstream <- mcols(
        tr_upstream_track[
          GenomicRanges::start(tr_upstream_track@range) ==
            fusion@gene_upstream@breakpoint
        ]@range
      )$exon
    } else {
      exon_breakpoint_upstream <- mcols(
        tr_upstream_track[
          GenomicRanges::end(tr_upstream_track@range) ==
            fusion@gene_upstream@breakpoint
        ]@range
      )$exon
    }
    if (gene_downstream_at_minus_strand) {
      exon_breakpoint_downstream <- mcols(
        tr_downstream_track[
          GenomicRanges::end(tr_downstream_track@range) ==
            fusion@gene_downstream@breakpoint
        ]@range
      )$exon
    } else {
      exon_breakpoint_downstream <- mcols(
        tr_downstream_track[
          GenomicRanges::start(tr_downstream_track@range) ==
            fusion@gene_downstream@breakpoint
        ]@range
      )$exon
    }

    # Get coordinates for the exons in the plot
    exon_coordinates1 <- Gviz::coords(res1$GeneRegionTrack)
    exon_coordinates2 <- Gviz::coords(res2$GeneRegionTrack)

    # Since we know which exon in geneA meets the breakpoint, its easy to fetch
    # the coordinates for the breakpoint exon in geneA:
    exon_upstream <- exon_coordinates1[
      exon_breakpoint_upstream,
      ,
      drop = FALSE
    ]
    # Do the same for geneB:
    exon_downstream <- exon_coordinates2[
      exon_breakpoint_downstream,
      ,
      drop = FALSE
    ]

    # The coordinates we're really after is the top right corner of exonA (or
    # top left corner if on the minus strand), and the top left (top right if on
    # minus strand) corner of exon B
    if (gene_upstream_at_minus_strand) {
      # top left corner:
      x_pos_gene_upstream <- min(exon_upstream[, 1])
      y_pos_gene_upstream <- min(exon_upstream[, 2])
    } else {
      # top right corner:
      x_pos_gene_upstream <- min(exon_upstream[, 3])
      y_pos_gene_upstream <- min(exon_upstream[, 2])
    }
    if (gene_downstream_at_minus_strand) {
      # top right corner:
      x_pos_gene_downstream <- min(exon_downstream[, 3])
      y_pos_gene_downstream <- min(exon_downstream[, 2])
    } else {
      # top left corner:
      x_pos_gene_downstream <- min(exon_downstream[, 1])
      y_pos_gene_downstream <- min(exon_downstream[, 2])
    }

    # This will keep the coordinates even though the dev.size change:
    grid::pushViewport(
      grid::viewport(
        xscale = c(
          0,
          grDevices::dev.size(units = "px")[1]
        ), # x goes from 0 to device width
        yscale = c(
          grDevices::dev.size(units = "px")[2],
          0
        ) # y goes from device height to 0
      )
    )

    # Calculate the bezier curve width. The following will have the effect of
    # setting the line width to 10 if there are 50 or more reads that support
    # the fusion. If there are between 1 and 50 reads supporting the fusion,
    # then the curve width will be scaled accordingly.
    supporting_reads <- fusion@spanning_reads_count + fusion@split_reads_count
    supporting_reads_text <- paste0(
      supporting_reads,
      " (",
      fusion@split_reads_count,
      ",",
      fusion@spanning_reads_count,
      ")")
    curve_width <- .scale_list_to_interval(
      c(supporting_reads,
        1,
        20),
      1,
      5)
    curve_width <- curve_width[1]

    # How far above the transcript should the bezier curve control points be?
    bezier_control_point_offset <- 18

    # The different transcripts will be plotted at different heights. Choose the
    # y position closest to the top of the plot. That is the lowest y value of
    # the two, since the y axis is flipped
    y_pos <- min(y_pos_gene_upstream, y_pos_gene_downstream)

    # Draw the bezier curve between the transcripts
    grid::grid.bezier(
      rep(
        c(
          x_pos_gene_upstream,
          x_pos_gene_downstream
        ),
        each = 2
      ), # x positions
      c(y_pos,
        y_pos - bezier_control_point_offset,
        y_pos - bezier_control_point_offset,
        y_pos
      ), # y positions
      default.units = "native",
      gp = grid::gpar(
        col = "red",
        lwd = curve_width))

    # Where should the number of supporting reads be plotted? Calculate position
    number_offset <- 7
    number_position_x <-
      x_pos_gene_upstream + (x_pos_gene_downstream - x_pos_gene_upstream) / 2
    number_position_y <-
      y_pos - bezier_control_point_offset - number_offset

    grid::grid.text(
      supporting_reads_text,
      x = number_position_x / grDevices::dev.size(units = "px")[1],
      y = 1 - number_position_y / grDevices::dev.size(units = "px")[2],
      vp = grid::viewport(
        xscale = c(0, grDevices::dev.size(units = "px")[1]),
        yscale = c(grDevices::dev.size(units = "px")[2], 0)),
      gp = grid::gpar(
        fontsize = 15))

  }
}

#' @name plot_fusion
#'
#' @export
plot_fusion_together <- function(
  fusion,
  edb,
  bamfile = NULL,
  which_transcripts = "exonBoundary",
  ylim = c(0, 1000),
  non_ucsc = TRUE,
  reduce_transcripts = FALSE,
  bedgraphfile = NULL) {

  .validate_plot_fusion_params(
    fusion,
    edb,
    bamfile,
    which_transcripts,
    ylim,
    non_ucsc,
    reduce_transcripts,
    bedgraphfile
  )
  fusion <- .get_transcripts_if_not_there(fusion, edb)

  # Create tracks

  # Select transcripts
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
    # Add symbol metadata column to make Gviz group correctly
    mcols(reduced_transcripts_upstream)$symbol <- fusion@gene_upstream@name
    mcols(reduced_transcripts_downstream)$symbol <- fusion@gene_downstream@name
    # Create the track
    tr_track_both <- Gviz::GeneRegionTrack(
      rbind(
        data.frame(reduced_transcripts_upstream),
        data.frame(reduced_transcripts_downstream)
      ),
      chromosome = fusion@gene_upstream@chromosome)
  } else {
    # Create the track
    tr_track_both <- Gviz::GeneRegionTrack(
      rbind(
        data.frame(transcripts_upstream),
        data.frame(transcripts_downstream)
      ),
      chromosome = fusion@gene_upstream@chromosome)
  }
  # Set display parameters
  transcript_display_params <- list(
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
  Gviz::displayPars(tr_track_both) <- transcript_display_params

  # Create axis track
  axis_track <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axis_track) <- list(
    add53 = TRUE,
    add35 = TRUE,
    col = "black", # line color
    fontcolor = "black",
    fontsize = 15
  )

  # Create ideogram tracks

  # Read cytoband information depending on genome version
  if (fusion@genome_version == "hg19") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.HG19.Human.CytoBandIdeogram.txt",
      package = "chimeraviz")
  } else if (fusion@genome_version == "hg38") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.HG38.Human.CytoBandIdeogram.txt",
      package = "chimeraviz")
  } else if (fusion@genome_version == "mm10") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.MM10.Mus.musculus.CytoBandIdeogram.txt",
      package = "chimeraviz")
  } else {
    stop("Unsupported genome version")
  }
  cytoband <- utils::read.table(cytoband_file)
  # Set names to what Gviz expects
  names(cytoband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

  ideogram_track <- Gviz::IdeogramTrack(
    genome = fusion@genome_version,
    chromosome = fusion@gene_upstream@chromosome,
    bands = cytoband)
  # Set display paramters
  ideogram_display_params <- list(
    showBandId = TRUE,
    showId = TRUE,
    fontcolor = "black",
    cex = 1.5)
  Gviz::displayPars(ideogram_track) <- ideogram_display_params

  # Create alignment track
  if (!is.null(bamfile)) {
    # We're getting coverage data from a bam file
    if (non_ucsc) {
      # If the bam file has non-ucsc chromosome names, i.e. "1" instead of
      # "chr1", then we need to use a custom import function that adds "chr" to
      # the chromosome names, to keep Gviz happy. Gviz strongly prefers having
      # the "chr" prefix.
      al_track <- Gviz::AlignmentsTrack(
        bamfile,
        isPaired = TRUE,
        genome = fusion@genome_version,
        name = "RNA coverage",
        ylim = ylim,
        importFunction = import_function_non_ucsc)
    } else {
      al_track <- Gviz::AlignmentsTrack(
        bamfile,
        isPaired = TRUE,
        genome = fusion@genome_version,
        name = "RNA coverage",
        ylim = ylim)
    }
    # Set display paramters
    Gviz::displayPars(al_track) <- list(
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
  } else {
    # We're getting coverage data from a bedGraph file
    al_track <- DataTrack(
      range = bedgraphfile,
      ylim = ylim,
      genome = "hg19",
      chromosome = "1",
      name = "Coverage",
      type = "h",
      col = "orange",
      fill = "orange")
    # Set display parameters
    Gviz::displayPars(al_track) <- list(
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


  # Display parameters for highlight tracks
  display_pars_highlight_tracks <- list(
    fill = "lightgrey",
    col = "lightgrey",
    showTitle = FALSE
  )

  # Helper variables that decide which way to plot the transcript later
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
  if (reduce_transcripts) {
    gene_upstream_indexes <-
      Gviz::symbol(tr_track_both) == fusion@gene_upstream@name
    gene_downstream_indexes <-
      Gviz::symbol(tr_track_both) == fusion@gene_downstream@name
  } else {
    gene_upstream_indexes <-
      Gviz::symbol(tr_track_both) %in% names(transcripts_upstream)
    gene_downstream_indexes <-
      Gviz::symbol(tr_track_both) %in% names(transcripts_downstream)
  }

  # Set all colors to exluded
  Gviz::feature(tr_track_both)[
    Gviz::feature(tr_track_both) == "5utr" & gene_upstream_indexes
  ] <- "5utr_excludedA"
  Gviz::feature(tr_track_both)[
    Gviz::feature(tr_track_both) == "protein_coding" & gene_upstream_indexes
  ] <- "protein_coding_excludedA"
  Gviz::feature(tr_track_both)[
    Gviz::feature(tr_track_both) == "3utr" & gene_upstream_indexes
  ] <- "3utr_excludedA"

  Gviz::feature(tr_track_both)[
    Gviz::feature(tr_track_both) == "5utr" & gene_downstream_indexes
  ] <- "5utr_excludedB"
  Gviz::feature(tr_track_both)[
    Gviz::feature(tr_track_both) == "protein_coding" & gene_downstream_indexes
  ] <- "protein_coding_excludedB"
  Gviz::feature(tr_track_both)[
    Gviz::feature(tr_track_both) == "3utr" & gene_downstream_indexes
  ] <- "3utr_excludedB"

  if (gene_upstream_at_minus_strand) {
    message(paste("As ",
                  fusion@gene_upstream@name,
                  " is on the minus strand, the plot for this gene will be ",
                  "reversed",
                  sep = ""))
    highlight_start_upstream <- fusion@gene_upstream@breakpoint
    highlight_end_upstream <-
      max(
        start(
          tr_track_both[gene_upstream_indexes]
        ) + width(tr_track_both[gene_upstream_indexes]) - 1)

    # Color the exons included/not included in the fusion
    Gviz::feature(tr_track_both)[
      start(tr_track_both) <= highlight_end_upstream &
      end(tr_track_both) >= highlight_start_upstream &
      Gviz::feature(tr_track_both) == "5utr_excludedA" &
      gene_upstream_indexes
    ] <- "5utr_includedA"

    Gviz::feature(tr_track_both)[
      start(tr_track_both) <= highlight_end_upstream &
      end(tr_track_both) >= highlight_start_upstream &
      Gviz::feature(tr_track_both) == "protein_coding_excludedA" &
      gene_upstream_indexes
    ] <- "protein_coding_includedA"

    Gviz::feature(tr_track_both)[
      start(tr_track_both) <= highlight_end_upstream &
      end(tr_track_both) >= highlight_start_upstream &
      Gviz::feature(tr_track_both) == "3utr_excludedA" &
      gene_upstream_indexes
    ] <- "3utr_includedA"
  } else {
    highlight_start_upstream <- min(start(tr_track_both[gene_upstream_indexes]))
    highlight_end_upstream <- fusion@gene_upstream@breakpoint

    # Color the exons included/not included in the fusion
    Gviz::feature(tr_track_both)[
      start(tr_track_both) >= highlight_start_upstream &
      end(tr_track_both) <= highlight_end_upstream &
      Gviz::feature(tr_track_both) == "5utr_excludedA" &
      gene_upstream_indexes
    ] <- "5utr_includedA"

    Gviz::feature(tr_track_both)[
      start(tr_track_both) >= highlight_start_upstream &
      end(tr_track_both) <= highlight_end_upstream &
      Gviz::feature(tr_track_both) == "protein_coding_excludedA" &
      gene_upstream_indexes
    ] <- "protein_coding_includedA"

    Gviz::feature(tr_track_both)[
      start(tr_track_both) >= highlight_start_upstream &
      end(tr_track_both) <= highlight_end_upstream &
      Gviz::feature(tr_track_both) == "3utr_excludedA" &
      gene_upstream_indexes
    ] <- "3utr_includedA"
  }

  if (gene_downstream_at_minus_strand) {
    message(paste("As ",
                  fusion@gene_downstream@name,
                  " is on the minus strand, the plot for this gene will be ",
                  "reversed",
                  sep = ""))
    highlight_start_downstream <-
      min(start(tr_track_both[gene_downstream_indexes]))
    highlight_end_downstream <- fusion@gene_downstream@breakpoint

    # Color the exons included/not included in the fusion
    Gviz::feature(tr_track_both)[
      start(tr_track_both) <= highlight_end_downstream &
      end(tr_track_both) >= highlight_start_downstream &
      Gviz::feature(tr_track_both) == "5utr_excludedB" &
      gene_downstream_indexes
    ] <- "5utr_includedB"

    Gviz::feature(tr_track_both)[
      start(tr_track_both) <= highlight_end_downstream &
      end(tr_track_both) >= highlight_start_downstream &
      Gviz::feature(tr_track_both) == "protein_coding_excludedB" &
      gene_downstream_indexes
    ] <- "protein_coding_includedB"

    Gviz::feature(tr_track_both)[
      start(tr_track_both) <= highlight_end_downstream &
      end(tr_track_both) >= highlight_start_downstream &
      Gviz::feature(tr_track_both) == "3utr_excludedB" &
      gene_downstream_indexes
    ] <- "3utr_includedB"
  } else {
    highlight_start_downstream <- fusion@gene_downstream@breakpoint
    highlight_end_downstream <-
      max(
        start(
          tr_track_both[gene_downstream_indexes]
        ) + width(tr_track_both[gene_downstream_indexes]) - 1)

    # Color the exons included/not included in the fusion
    Gviz::feature(tr_track_both)[
      start(tr_track_both) >= highlight_start_downstream &
      end(tr_track_both) <= highlight_end_downstream &
      Gviz::feature(tr_track_both) == "5utr_excludedB" &
      gene_downstream_indexes
    ] <- "5utr_includedB"

    Gviz::feature(tr_track_both)[
      start(tr_track_both) >= highlight_start_downstream &
      end(tr_track_both) <= highlight_end_downstream &
      Gviz::feature(tr_track_both) == "protein_coding_excludedB" &
      gene_downstream_indexes
    ] <- "protein_coding_includedB"

    Gviz::feature(tr_track_both)[
      start(tr_track_both) >= highlight_start_downstream &
      end(tr_track_both) <= highlight_end_downstream &
      Gviz::feature(tr_track_both) == "3utr_excludedB" &
      gene_downstream_indexes
    ] <- "3utr_includedB"
  }

  # Update transcript track with colors
  Gviz::displayPars(tr_track_both) <- interestcolor

  # Highlight transcript track

  gr_track_highlight <- Gviz::HighlightTrack(
    trackList = tr_track_both,
    start = c(highlight_start_upstream, highlight_start_downstream),
    end = c(highlight_end_upstream, highlight_end_downstream),
    chromosome = fusion@gene_upstream@chromosome)
  Gviz::displayPars(gr_track_highlight) <- display_pars_highlight_tracks

  # Highlight coverage track

  al_track_highlight <- Gviz::HighlightTrack(
    trackList = al_track,
    start = c(highlight_start_upstream, highlight_start_downstream),
    end = c(highlight_end_upstream, highlight_end_downstream),
    chromosome = fusion@gene_upstream@chromosome)
  Gviz::displayPars(al_track_highlight) <- display_pars_highlight_tracks

  # Do some plotting

  # Open a new page to plot
  grid::grid.newpage()

  # How tall should the transcripts row be? This depends on how many transcripts
  # we want to show. If reduce_transcripts=TRUE, then we only want to display
  # one transcript per gene.
  if (reduce_transcripts) {
    transcripts_row_height <- 1
  } else {
    # It depends on max(amountOfTranscriptsA, amountOfTranscriptsB)
    transcripts_row_height <- max(
      length(transcripts_upstream),
      length(transcripts_downstream))
    # Set a max value
    if (transcripts_row_height > 10) {
      transcripts_row_height <- 10
    }
  }

  # If the range for the plot is too short, then GenomeAxisTrack will show
  # ticks too close to each other. This doesn't look nice. A way to solve it is
  # to extend the plot range if the transcript is too short:
  plot_range <- max(end(tr_track_both)) - min(start(tr_track_both))
  if (plot_range < 20000) {
    offset <- 2500
  } else {
    offset <- 0
  }

  # Plot transcripts and coverage with highlight
  res <- Gviz::plotTracks(
    # without this gviz create cluster_X entries in the GeneRegionTrack
    collapse = FALSE,
    c(ideogram_track, gr_track_highlight, al_track_highlight, axis_track),
    # 10k added so that we can see all transcript names
    from = min(start(tr_track_both)) - offset - 10000,
    # 10k added so that we can see all transcript names
    to = max(end(tr_track_both)) + offset + 10000,
    sizes = c(1, transcripts_row_height, 2, 1),
    add = TRUE,
    background.title = "transparent",
    chromosome = rep(fusion@gene_upstream@chromosome, 4),
    # Plot reverse if gene is at minus strand
    reverseStrand = gene_upstream_at_minus_strand)

  # Draw bezier curve if we're plotting exon boundary transcripts

  upstream_boundary_exons <-
    fusion@gene_upstream@transcripts[
      mcols(fusion@gene_upstream@transcripts)$transcript_category ==
        "exonBoundary"
    ]
  downstream_boundary_exons <-
    fusion@gene_downstream@transcripts[
      mcols(fusion@gene_downstream@transcripts)$transcript_category ==
        "exonBoundary"
    ]
  if (length(upstream_boundary_exons) > 0 &
      length(downstream_boundary_exons) > 0) {

    # We need to figure out which exons in trTrackBoth matches the
    # fusion breakpoints fusion@gene_upstream@breakpoint and
    # fusion@gene_downstream@breakpoint:
    if (gene_upstream_at_minus_strand) {
      exon_breakpoint_upstream <-
        mcols(
          tr_track_both[
            GenomicRanges::start(tr_track_both@range) ==
              fusion@gene_upstream@breakpoint
            ]@range
        )$exon
    } else {
      exon_breakpoint_upstream <-
        mcols(
          tr_track_both[
            GenomicRanges::end(tr_track_both@range) ==
              fusion@gene_upstream@breakpoint
            ]@range
        )$exon
    }
    if (gene_downstream_at_minus_strand) {
      exon_breakpoint_downstream <-
        mcols(
          tr_track_both[
            GenomicRanges::end(tr_track_both@range) ==
              fusion@gene_downstream@breakpoint
          ]@range
        )$exon
    } else {
      exon_breakpoint_downstream <-
        mcols(
          tr_track_both[
            GenomicRanges::start(tr_track_both@range) ==
              fusion@gene_downstream@breakpoint
          ]@range
        )$exon
    }

    # Get coordinates for the exons in the plot
    exo_coordinates <- Gviz::coords(res$GeneRegionTrack)

    # Since we know which exon in geneA meets the breakpoint, its easy to fetch
    # the coordinates for the breakpoint exon in geneA:
    exon_upstream <- exo_coordinates[
      exon_breakpoint_upstream,
      ,
      drop = FALSE
    ]
    # Do the same for geneB:
    exon_downstream <- exo_coordinates[
      exon_breakpoint_downstream,
      ,
      drop = FALSE
    ]

    # The coordinates we're really after is the top right corner of exonA (or
    # top left corner if on the minus strand), and the top left (top right if on
    # minus strand) corner of exon B
    if (gene_upstream_at_minus_strand) {
      # top left corner:
      x_pos_gene_upstream <- min(exon_upstream[, 1])
      y_pos_gene_upstream <- min(exon_upstream[, 2])
    } else {
      # top right corner:
      x_pos_gene_upstream <- min(exon_upstream[, 3])
      y_pos_gene_upstream <- min(exon_upstream[, 2])
    }
    if (gene_downstream_at_minus_strand) {
      # top right corner:
      x_pos_gene_downstream <- min(exon_downstream[, 3])
      y_pos_gene_downstream <- min(exon_downstream[, 2])
    } else {
      # top left corner:
      x_pos_gene_downstream <- min(exon_downstream[, 1])
      y_pos_gene_downstream <- min(exon_downstream[, 2])
    }

    # This will keep the coordinates even though the dev.size change:
    grid::pushViewport(
      grid::viewport(
        xscale = c(
          0,
          grDevices::dev.size(units = "px")[1]
        ), # x goes from 0 to device width
        yscale = c(
          grDevices::dev.size(units = "px")[2],
          0
        ) # y goes from device height to 0
      )
    )

    # # Draw some point just to check that we have the right coordinates:
    # grid.points(topCornerExonA[1], topCornerExonA[2]) # Exclude Linting
    # grid.points(topCornerExonB[1], topCornerExonB[2]) # Exclude Linting

    # Calculate the bezier curve width. The following will have the effect of
    # setting the line width to 10 if there are 50 or more reads that support
    # the fusion. If there are between 1 and 50 reads supporting the fusion,
    # then the curve width will be scaled accordingly.
    supporting_reads <- fusion@spanning_reads_count + fusion@split_reads_count
    supporting_reads_text <- paste0(
      supporting_reads,
      " (",
      fusion@split_reads_count,
      ",",
      fusion@spanning_reads_count,
      ")")
    curve_width <- .scale_list_to_interval(
      c(supporting_reads,
        1,
        20),
      1,
      5)
    curve_width <- curve_width[1]

    # How far above the transcript should the bezier curve control points be?
    bezier_control_point_offset <- 18

    # The different transcripts will be plotted at different heights. Choose the
    # y position closest to the top of the plot. That is the lowest y value of
    # the two, since the y axis is flipped
    y_pos <- min(y_pos_gene_upstream, y_pos_gene_downstream)

    # Draw the bezier curve between the transcripts
    grid::grid.bezier(
      rep(
        c(
          x_pos_gene_upstream,
          x_pos_gene_downstream
        ),
        each = 2
      ), # x positions
      c(
        y_pos,
        y_pos - bezier_control_point_offset,
        y_pos - bezier_control_point_offset,
        y_pos
      ), # y positions
      default.units = "native",
      gp = grid::gpar(
        col = "red",
        lwd = curve_width))

    # Where should the number of supporting reads be plotted? Calculate position
    number_offset <- 7
    number_position_x <-
      x_pos_gene_upstream + (x_pos_gene_downstream - x_pos_gene_upstream) / 2
    number_position_y <-
      y_pos - bezier_control_point_offset - number_offset

    grid::grid.text(
      supporting_reads_text,
      x = number_position_x / grDevices::dev.size(units = "px")[1],
      y = 1 - number_position_y / grDevices::dev.size(units = "px")[2],
      vp = grid::viewport(
        xscale = c(0, grDevices::dev.size(units = "px")[1]),
        yscale = c(grDevices::dev.size(units = "px")[2], 0)),
      gp = grid::gpar(
        fontsize = 15))

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
import_function_non_ucsc <- function (file, selection) {

  ind_names <- c(sub("\\.bam$", ".bai", file), paste(file,
                                                    "bai", sep = "."))
  index <- NULL
  for (i in ind_names) {
    if (file.exists(i)) {
      index <- i
      break
    }
  }
  if (is.null(index))
    stop(
      "Unable to find index for BAM file '",
      file,
      "'. You can build an index using the following command:\n\t",
      "library(Rsamtools)\n\tindexBam(\"",
      file,
      "\")"
    )
  paired_end <- parent.env(environment())[["._isPaired"]]
  if (is.null(paired_end))
    paired_end <- TRUE
  bf <- Rsamtools::BamFile(file, index = index, asMates = paired_end)
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
  reads <- if (
    as.character(seqnames(selection)[1]) %in%
    names(Rsamtools::scanBamHeader(bf)$targets)
  ) {
    Rsamtools::scanBam(bf, param = param)[[1]]
  } else {
    list()
  }
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
  return(
    GRanges(
      # Return chromosome name prefixed with "chr"
      seqnames = if (is.null(reads$rname)) {
        character()
      } else {
        paste("chr", reads$rname, sep = "")
      },
      strand = if (is.null(reads$strand)) {
        character()
      } else {
        reads$strand
      },
      ranges = IRanges(start = reads$pos, width = reads$qwidth),
      id = reads$qname, cigar = reads$cigar, mapq = reads$mapq,
      flag = reads$flag, md = md, seq = ans, isize = reads$isize,
      groupid = if (paired_end) {
        reads$groupid
      } else {
        seq_along(reads$pos)
      },
      status = if (paired_end) {
        reads$mate_status
      } else {
        rep(
          factor(
            "unmated",
            levels = c("mated", "ambiguous", "unmated")
          ),
          length(reads$pos)
        )
      }
    )
  )
}

.validate_plot_fusion_params <- function(
  fusion,
  edb = NULL,
  bamfile,
  which_transcripts = "exonBoundary",
  ylim = c(0, 1000),
  non_ucsc = TRUE,
  reduce_transcripts = FALSE,
  bedgraphfile
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()

  # Check parameters
  argument_checker <- .is_fusion_valid(argument_checker, fusion)
  argument_checker <- .is_edb_valid(argument_checker, edb, fusion)
  argument_checker <- .is_either_bamfile_or_bedgraphfile_valid(
    argument_checker,
    bamfile,
    bedgraphfile
  )
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
