#' Graph plot of possible fusion transcripts.
#'
#' This function takes a fusion object and a TranscriptDb object and plots a
#' graph showing the possible fusion transcripts.
#'
#' @param fusion The Fusion object to plot.
#' @param edb The edb object that will be used to fetch data.
#' @param which_transcripts This character vector decides which transcripts are
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
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
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
#' plot_fusion_transcripts_graph(
#'   fusion = fusion,
#'   edb = edb)
#' # Close device
#' dev.off()
#'
#' @export
plot_fusion_transcripts_graph <- function(
  fusion,
  edb = NULL,
  which_transcripts = "exonBoundary",
  rankdir = "TB") {

  .validate_plot_fusion_transcripts_graph_params(
    fusion,
    edb,
    which_transcripts,
    rankdir
  )
  fusion <- .get_transcripts_if_not_there(fusion, edb)

  # Select which transcripts to use
  transcripts_upstream <-
    select_transcript(fusion@gene_upstream, which_transcripts)
  transcripts_downstream <-
    select_transcript(fusion@gene_downstream, which_transcripts)

  # Unlist
  transcripts_upstream <- unlist(transcripts_upstream)
  transcripts_downstream <- unlist(transcripts_downstream)

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

  if (fusion@gene_upstream@strand == "+") {
    fusion_exons_upstream <-
      transcripts_upstream[
        start(ranges(transcripts_upstream)) < fusion@gene_upstream@breakpoint &
          end(ranges(transcripts_upstream)) <= fusion@gene_upstream@breakpoint
      ]
  } else {
    fusion_exons_upstream <-
      transcripts_upstream[
        start(ranges(transcripts_upstream)) >= fusion@gene_upstream@breakpoint &
          end(ranges(transcripts_upstream)) > fusion@gene_upstream@breakpoint
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

  # If either fusionExonsA or fusionExonsB is empty, then the breakpoint does
  # not occur inside any of the transcripts.
  if (length(fusion_exons_upstream) == 0) {
    stop(paste(
      "None of the transcripts given for gene A has the fusion breakpoint ",
      "within them. This plot cannot be created with the given transcripts."))
  }
  if (length(fusion_exons_downstream) == 0) {
    stop(paste(
      "None of the transcripts given for gene B has the fusion breakpoint ",
      "within them. This plot cannot be created with the given transcripts."))
  }

  # For each transcript in fusionExonsA, match the exons in those transcripts..

  transcript_names_upstream <-
    unique(unlist(mcols(fusion_exons_upstream)$transcript))
  transcript_names_downstream <-
    unique(unlist(mcols(fusion_exons_downstream)$transcript))
  fusion_transcripts <- GRangesList()
  for (i in seq_along(transcript_names_upstream)) {
    tr_name_upstream <- transcript_names_upstream[[i]]
    for (j in seq_along(transcript_names_downstream)) {
      tr_name_downstream <- transcript_names_downstream[[j]]

      exons_upstream <- fusion_exons_upstream[
        mcols(fusion_exons_upstream)$transcript == tr_name_upstream
      ]
      exons_downstream <- fusion_exons_downstream[
        mcols(fusion_exons_downstream)$transcript == tr_name_downstream
      ]

      # Reverse order of exons if on minus strand
      if (fusion@gene_upstream@strand == "-") {
        exons_upstream <- rev(exons_upstream)
      }
      if (fusion@gene_downstream@strand == "-") {
        exons_downstream <- rev(exons_downstream)
      }

      new <- append(
        exons_upstream,
        exons_downstream
      )
      mcols(new)$transcript <- paste("fusion", i, j, sep = "")
      fusion_transcripts <- append(
        fusion_transcripts,
        GRangesList(new)
      )
    }
  }

  fusion_transcripts <- unlist(fusion_transcripts)

  exon_names_upstream <-
    unique(
      mcols(
        fusion_transcripts[
          mcols(fusion_transcripts)$gene_id == fusion@gene_upstream@ensembl_id
        ]
      )$exon_id
    )
  exon_names_downstream <-
    unique(
      mcols(
        fusion_transcripts[
          mcols(fusion_transcripts)$gene_id == fusion@gene_downstream@ensembl_id
        ]
      )$exon_id
    )

  # These lists will hold edge information
  from <- c()
  to <- c()

  transcript_names <- sort(unique(mcols(fusion_transcripts)$transcript))
  for (i in seq_along(transcript_names)) {
    exons_in_transcript <-
      fusion_transcripts[
        mcols(fusion_transcripts)$transcript == transcript_names[i]
      ]
    for (j in 2:length(exons_in_transcript)) {
      # Add normal edge
      from <- c(from, mcols(exons_in_transcript[j - 1])$exon_id)
      to <- c(to, mcols(exons_in_transcript[j])$exon_id)
    }
  }

  # Use the edge information we've gathered and create a graph
  ft <- cbind(from, to)
  # Remove edges that go to the same node. Do this because the way
  # get_transcripts_ensembl_db is implemented it will split exons in which the
  # cds starts or ends. The two new exons from such a split will have the same
  # id, causing these edges to appear.
  ft <- ft[ft[, 1] != ft[, 2], ]
  ft[ft[, 1] == ft[, 2], ]

  r_eg <- graph::ftM2graphNEL(unique(ft))

  # Color exons
  cols <- RColorBrewer::brewer.pal(4, "Paired")
  interestcolor <- list(
    "excludedA" = cols[1],
    "includedA" = cols[2],
    "excludedB" = cols[3],
    "includedB" = cols[4])

  # Create a list of exon names with colors
  list_upstream <- rep(interestcolor[2], length(exon_names_upstream))
  names(list_upstream) <- exon_names_upstream
  list_downstream <- rep(interestcolor[4], length(exon_names_downstream))
  names(list_downstream) <- exon_names_downstream

  n_attrs <- list()
  e_attrs <- list()
  n_attrs$fillcolor <- c(
    "start" = "gray",
    list_upstream,
    list_downstream,
    "stop" = "gray")

  Rgraphviz::plot(
    r_eg,
    attrs = list(
      node = list(
        shape = "box",
        fontsize = 42),
      edge = list(),
      graph = list(rankdir = rankdir)), # left-to-right: LR, top-to-bottom: TB
    nodeAttrs = n_attrs,
    edgeAttrs = e_attrs)

}

.validate_plot_fusion_transcripts_graph_params <- function(
  fusion,
  edb,
  which_transcripts,
  rankdir
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()

  # Check parameters
  argument_checker <- .is_fusion_valid(argument_checker, fusion)
  argument_checker <- .is_edb_valid(argument_checker, edb, fusion)
  argument_checker <- .is_which_transcripts_valid(
    argument_checker,
    which_transcripts,
    fusion)

  # The "withinExon" option for which_transcripts will not work for the fusion
  # transcripts graph plot. Check if that was the wanted option
  if (which_transcripts[[1]] == "withinExon") {
    ArgumentCheck::addError(
      msg = paste0("The \"withinExon\" option for which_transcripts will not ",
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
