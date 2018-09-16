#' Plot a specific fusion transcript with protein domain annotations
#'
#' This function takes a fusion object, an ensembldb object, a bedfile with
#' protein domain data and two specific transcript ids. The function plots the
#' specific fusion transcript along with annotations of protein domains. If a
#' bamfile is specified, the fusion transcript will be plotted with coverage
#' information.
#'
#' Note that the transcript database used (the edb object) must have the same
#' seqnames as any bamfile used. Otherwise the coverage data will be wrong.
#'
#' @param fusion The Fusion object to plot.
#' @param edb The edb object that will be used to fetch data.
#' @param bamfile The bamfile with RNA-seq data.
#' @param bedfile The bedfile with protein domain data.
#' @param gene_upstream_transcript The transcript id for the upstream gene.
#' @param gene_downstream_transcript The transcript id for the downstream gene.
#' @param plot_downstream_protein_domains_if_fusion_is_out_of_frame Setting this to true
#' makes the function plot protein domains in the downstream gene even though
#' the fusion is out of frame.
#'
#' @return Creates a fusion transcript plot with annotations of protein domains.
#'
#' @examples
#' # Load data and example fusion event
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Select transcripts
#' gene_upstream_transcript <- "ENST00000434290"
#' gene_downstream_transcript <- "ENST00000370031"
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
#' # bedfile with protein domains for the transcripts in this example
#' bedfile <- system.file(
#'   "extdata",
#'   "protein_domains_5267.bed",
#'   package="chimeraviz")
#' # Temporary file to store the plot
#' pngFilename <- tempfile(
#'   pattern = "fusionPlot",
#'   fileext = ".png",
#'   tmpdir = tempdir())
#' # Open device
#' png(pngFilename, width = 500, height = 500)
#' # Plot!
#' plot_fusion_transcript_with_protein_domain(
#'   fusion = fusion,
#'   edb = edb,
#'   bamfile = bamfile5267,
#'   bedfile = bedfile,
#'   gene_upstream_transcript = gene_upstream_transcript,
#'   gene_downstream_transcript = gene_downstream_transcript,
#'   plot_downstream_protein_domains_if_fusion_is_out_of_frame = TRUE)
#' # Close device
#' dev.off()
#'
#' @export
#'
plot_fusion_transcript_with_protein_domain <- function(
  fusion,
  edb = NULL,
  bamfile = NULL,
  bedfile = NULL,
  gene_upstream_transcript = "",
  gene_downstream_transcript = "",
  plot_downstream_protein_domains_if_fusion_is_out_of_frame = FALSE) {

  # Validate that we're given a bedfile that exists
  if (class(bedfile) != "character") {
    stop("bedfile argument must be a character")
  }
  if (file.exists(bedfile) == FALSE) {
    stop("Invalid bedfile: File doesn't exist.")
  }
  # Validate geneAtranscript and geneBtranscript
  if (class(gene_upstream_transcript) != "character") {
    stop("gene_upstream_transcript argument must be a character")
  }
  if (class(gene_downstream_transcript) != "character") {
    stop("gene_downstream_transcript argument must be a character")
  }
  if (gene_upstream_transcript == "") {
    stop("gene_upstream_transcript argument must not be blank")
  }
  if (gene_downstream_transcript == "") {
    stop("gene_downstream_transcript argument must not be blank")
  }
  # Validate parameters and fetch transcripts from database
  .validate_plot_fusion_transcript_with_protein_domain_params(
    fusion,
    edb,
    bamfile,
    bedfile,
    gene_upstream_transcript,
    gene_downstream_transcript,
    plot_downstream_protein_domains_if_fusion_is_out_of_frame)
  fusion <- .get_transcripts_if_not_there(fusion, edb)
  # Validate that the transcript ids are found in the fusion
  if (
    length(
      fusion@gene_upstream@transcripts[
        names(fusion@gene_upstream@transcripts) == gene_upstream_transcript
      ]) == 0) {
    stop(
      paste0(
        gene_upstream_transcript,
        " was not found in among the transcripts for the upstream gene ",
        "partner."
      )
    )
  }
  if (
    length(
      fusion@gene_downstream@transcripts[
        names(fusion@gene_downstream@transcripts) == gene_downstream_transcript
      ]) == 0) {
    stop(
      paste0(
        gene_downstream_transcript,
        " was not found in among the transcripts for the downstream gene ",
        "partner."
      )
    )
  }
  # Extract the transcripts
  transcript_upstream <- unlist(
    fusion@gene_upstream@transcripts[
      names(fusion@gene_upstream@transcripts) == gene_upstream_transcript
    ]
  )
  transcript_downstream <- unlist(
    fusion@gene_downstream@transcripts[
      names(fusion@gene_downstream@transcripts) == gene_downstream_transcript
    ]
  )
  # Check that there's at least one transcript that has the fusion breakpoint
  # within the transcript.
  .check_that_breakpoints_are_within_transcripts(
    fusion,
    transcript_upstream,
    transcript_downstream
  )
  # Check that the transcripts have a breakpoint exon
  .check_that_transcripts_have_breakpoint_exons(
    fusion,
    transcript_upstream,
    transcript_downstream
  )
  # Add transcript metadata column to make Gviz group correctly
  mcols(transcript_upstream)$transcript <- fusion@gene_upstream@name
  mcols(transcript_downstream)$transcript <- fusion@gene_downstream@name
  # Add symbol metadata column to make Gviz show ids
  mcols(transcript_upstream)$symbol <- fusion@gene_upstream@name
  mcols(transcript_downstream)$symbol <- fusion@gene_downstream@name

  # Downshift transcripts (remove introns from the GRanges)
  transcript_upstream_downshifted <- down_shift(transcript_upstream)
  transcript_downstream_downshifted <- down_shift(transcript_downstream)

  # Figure out which exons from geneA and geneB are part of the fusion
  # transcript.
  #
  # Also, "rediscover" where the fusion breakpoints are in the downshifted
  # versions of the transcripts.
  #
  # 1: geneA is on the + strand
  #    - Exons part of the fusion transcript are those with start/end positions
  #      before breakpoint
  #    - Fusion breakpoint in the downshifted transcript is at the end
  #      coordinate of exon_id_at_fusion_breakpoint_a
  # 2: geneA is on the - strand
  #    - Exons part of the fusion transcript are those with start/end positions
  #      after breakpoint
  #    - Fusion breakpoint in the downshifted transcript is at the start
  #      coordinate of exon_id_at_fusion_breakpoint_a
  # 3: geneB is on the + strand
  #    - Exons part of the fusion transcript are those with start/end positions
  #      after breakpoint
  #    - Fusion breakpoint in the downshifted transcript is at the start
  #      coordinate of exon_id_at_fusion_breakpoint_b
  # 4: geneB is on the - strand
  #    - Exons part of the fusion transcript are those with start/end positions
  #      before breakpoint
  #    - Fusion breakpoint in the downshifted transcript is at the end
  #      coordinate of exon_id_at_fusion_breakpoint_b
  #
  if (fusion@gene_upstream@strand == "+") {

    exon_ids_part_of_fusion_from_transcript_a <- mcols(
      transcript_upstream[
        start(ranges(transcript_upstream)) < fusion@gene_upstream@breakpoint &
          end(ranges(transcript_upstream)) <= fusion@gene_upstream@breakpoint
      ]
    )$exon_id
    exons_part_of_fusion_a <- transcript_upstream[
      start(ranges(transcript_upstream)) < fusion@gene_upstream@breakpoint &
        end(ranges(transcript_upstream)) <= fusion@gene_upstream@breakpoint
    ]
    exon_id_at_fusion_breakpoint_a <- mcols(
      transcript_upstream[
        end(ranges(transcript_upstream)) == fusion@gene_upstream@breakpoint
      ]
    )$exon_id

    fusion_breakpoint_in_downshifted_transcript_a <-
      max(
        end(
          transcript_upstream_downshifted[
            mcols(transcript_upstream_downshifted)$exon_id ==
              exon_id_at_fusion_breakpoint_a
          ]
        )
      )
    exons_part_of_fusion_downshifted_a <-
      transcript_upstream_downshifted[
        start(ranges(transcript_upstream_downshifted)) <
          fusion_breakpoint_in_downshifted_transcript_a &
        end(ranges(transcript_upstream_downshifted)) <=
          fusion_breakpoint_in_downshifted_transcript_a
      ]

  } else {

    exon_ids_part_of_fusion_from_transcript_a <-
      mcols(
        transcript_upstream[
          start(ranges(transcript_upstream)) >=
            fusion@gene_upstream@breakpoint &
          end(ranges(transcript_upstream)) >
            fusion@gene_upstream@breakpoint
        ]
      )$exon_id
    exons_part_of_fusion_a <- transcript_upstream[
      start(ranges(transcript_upstream)) >=
        fusion@gene_upstream@breakpoint &
      end(ranges(transcript_upstream)) >
        fusion@gene_upstream@breakpoint
    ]
    exon_id_at_fusion_breakpoint_a <-
      mcols(
        transcript_upstream[
          start(ranges(transcript_upstream)) == fusion@gene_upstream@breakpoint
        ]
      )$exon_id

    fusion_breakpoint_in_downshifted_transcript_a <-
      min(
        start(
          transcript_upstream_downshifted[
            mcols(transcript_upstream_downshifted)$exon_id ==
              exon_id_at_fusion_breakpoint_a
          ]
        )
      )
    exons_part_of_fusion_downshifted_a <-
      transcript_upstream_downshifted[
        start(ranges(transcript_upstream_downshifted)) >=
          fusion_breakpoint_in_downshifted_transcript_a &
        end(ranges(transcript_upstream_downshifted)) >
          fusion_breakpoint_in_downshifted_transcript_a
      ]

  }
  if (fusion@gene_downstream@strand == "+") {

    exon_ids_part_of_fusion_from_transcript_b <- mcols(
      transcript_downstream[
        start(ranges(transcript_downstream)) >=
          fusion@gene_downstream@breakpoint &
        end(ranges(transcript_downstream)) >
          fusion@gene_downstream@breakpoint
      ]
    )$exon_id
    exons_part_of_fusion_b <- transcript_downstream[
      start(ranges(transcript_downstream)) >=
        fusion@gene_downstream@breakpoint &
      end(ranges(transcript_downstream)) >
        fusion@gene_downstream@breakpoint
    ]
    exon_id_at_fusion_breakpoint_b <- mcols(
      transcript_downstream[
        start(ranges(transcript_downstream)) ==
          fusion@gene_downstream@breakpoint
      ]
    )$exon_id

    fusion_breakpoint_in_downshifted_transcript_b <- min(
      start(
        transcript_downstream_downshifted[
          mcols(transcript_downstream_downshifted)$exon_id ==
            exon_id_at_fusion_breakpoint_b
        ]
      )
    )
    exons_part_of_fusion_downshifted_b <- transcript_downstream_downshifted[
      start(ranges(transcript_downstream_downshifted)) >=
        fusion_breakpoint_in_downshifted_transcript_b &
      end(ranges(transcript_downstream_downshifted)) >
        fusion_breakpoint_in_downshifted_transcript_b
    ]

  } else {

    exon_ids_part_of_fusion_from_transcript_b <- mcols(
      transcript_downstream[
        start(ranges(transcript_downstream)) <
          fusion@gene_downstream@breakpoint &
        end(ranges(transcript_downstream)) <=
          fusion@gene_downstream@breakpoint
      ]
    )$exon_id
    exons_part_of_fusion_b <- transcript_downstream[
      start(ranges(transcript_downstream)) <
        fusion@gene_downstream@breakpoint &
      end(ranges(transcript_downstream)) <=
        fusion@gene_downstream@breakpoint
    ]
    exon_id_at_fusion_breakpoint_b <- mcols(
      transcript_downstream[
        end(ranges(transcript_downstream)) == fusion@gene_downstream@breakpoint
      ]
    )$exon_id

    fusion_breakpoint_in_downshifted_transcript_b <- max(
      end(
        transcript_downstream_downshifted[
          mcols(transcript_downstream_downshifted)$exon_id ==
            exon_id_at_fusion_breakpoint_b
        ]
      )
    )
    exons_part_of_fusion_downshifted_b <- transcript_downstream_downshifted[
      start(ranges(transcript_downstream_downshifted)) <
        fusion_breakpoint_in_downshifted_transcript_b &
      end(ranges(transcript_downstream_downshifted)) <=
        fusion_breakpoint_in_downshifted_transcript_b
    ]

  }

  # Load protein domain data
  # Specify which columns from the bedfile we're going to read
  col_types_protein_data <- c(
    "Transcript_id" = "character",
    "Pfam_id" = "character",
    "Start" = "integer",
    "End" = "integer",
    "Domain_name_abbreviation" = "character",
    "Domain_name_full" = "character"
  )
  # Read the bedfile
  protein_data <- data.table::fread(
    input = bedfile,
    colClasses = col_types_protein_data,
    showProgress = FALSE
  )
  # Validate that the transcript ids are found in the bedfile
  if (
    nrow(
      protein_data[protein_data$Transcript_id == gene_upstream_transcript, ]
    ) == 0
  ) {
    stop(
      paste0(
        "Found no transcripts matching ",
        gene_upstream_transcript,
        " in the edb database."
      )
    )
  }
  if (
    nrow(
      protein_data[protein_data$Transcript_id == gene_downstream_transcript, ]
    ) == 0
  ) {
    stop(
      paste0(
        "Found no transcripts matching ",
        gene_downstream_transcript,
        " in the edb database."
      )
    )
  }
  # Extract the rows we need
  protein_data_gene_upstream <-
    protein_data[protein_data$Transcript_id == gene_upstream_transcript, ]
  protein_data_gene_downstream <-
    protein_data[protein_data$Transcript_id == gene_downstream_transcript, ]

  # For each protein domain we have (for each row in proteinDataGeneA and
  # proteinDataGeneb), decide which of the three categories the protein domain
  # is in: inside, outside, insideButBroken.
  #
  # The categories are explained as follows:
  # - inside:       Within the exons that are part of the fusion transcript.
  # - outside:      Not within the exons that are part of the fusion transcript
  # - insideButBroken: Partly within the exons that are part of the fusion
  #                 transcript (runs over one of the edges...)
  #
  # If a protein domain is inside, then the protein domain is retained in the
  # fusion transcript, and we can plot the annotation without problems.
  #
  # If a protein domain is outside, then it is not retained in the fusion
  # transcript and we can ignore it.
  #
  # If a protein domain is insideButBroken, then it is worth noting and we will
  # include it in the plot, but colored so that it's apparent that it is not
  # retained in the fusion transcript.
  #
  # Begin Exclude Linting
  # 1: geneA is on the + strand
  #    - The protein domain is "inside" if:
  #        protein_start_coordinate >=
  #          min(start(exons_part_of_fusion_downshifted_a)) &
  #        protein_end_coordinate <=
  #          max(end(exons_part_of_fusion_downshifted_a))
  #    - The protein domain is "insideButBroken" if:
  #        protein_start_coordinate >=
  #          min(start(exons_part_of_fusion_downshifted_a)) &
  #        protein_start_coordinate <
  #          max(end(exons_part_of_fusion_downshifted_a))
  #
  # 2: geneA is on the - strand
  #    - The protein domain is "inside" if:
  #        protein_start_coordinate <=
  #          max(end(exons_part_of_fusion_downshifted_a)) &
  #        protein_end_coordinate >=
  #          min(start(exons_part_of_fusion_downshifted_a))
  #    - The protein domain is "insideButBroken" if:
  #        protein_start_coordinate <=
  #          max(end(exons_part_of_fusion_downshifted_a)) &
  #        protein_start_coordinate >
  #          min(start(exons_part_of_fusion_downshifted_a))
  #
  # 3: geneB is on the + strand
  #    - The protein domain is "inside" if:
  #        protein_start_coordinate >=
  #          min(start(exons_part_of_fusion_downshifted_b)) &
  #        protein_end_coordinate <=
  #          max(end(exons_part_of_fusion_downshifted_b))
  #    - The protein domain is "insideButBroken" if:
  #        protein_start_coordinate >=
  #          min(start(exons_part_of_fusion_downshifted_b)) &
  #        protein_start_coordinate <
  #          max(end(exons_part_of_fusion_downshifted_b))
  #
  # 4: geneB is on the - strand
  #    - The protein domain is "inside" if:
  #        protein_start_coordinate <=
  #          max(end(exons_part_of_fusion_downshifted_b)) &
  #        protein_end_coordinate >=
  #          min(start(exons_part_of_fusion_downshifted_b))
  #    - The protein domain is "insideButBroken" if:
  #        protein_start_coordinate <=
  #          max(end(exons_part_of_fusion_downshifted_b)) &
  #        protein_start_coordinate >
  #          min(start(exons_part_of_fusion_downshifted_b))
  # End Exclude Linting

  # Go through all protein domains to classify them.
  # First extend the proteinDataGeneA and proteinDataGeneB data frames:
  protein_data_gene_upstream$protein_domain_location <- "none"
  protein_data_gene_upstream$plot_start <- "-"
  protein_data_gene_upstream$plot_end <- "-"
  protein_data_gene_downstream$protein_domain_location <- "none"
  protein_data_gene_downstream$plot_start <- "-"
  protein_data_gene_downstream$plot_end <- "-"
  # a)
  for (i in 1:nrow(protein_data_gene_upstream)) {
    if (
      protein_data_gene_upstream[i, ]$Start >=
        min(start(exons_part_of_fusion_downshifted_a)) &
      protein_data_gene_upstream[i, ]$End <=
        max(end(exons_part_of_fusion_downshifted_a))
    ) {
      # Is it "inside"?
      protein_data_gene_upstream[i, ]$protein_domain_location <- "inside"
      protein_data_gene_upstream[i, ]$plot_start <-
        protein_data_gene_upstream[i, ]$Start
      protein_data_gene_upstream[i, ]$plot_end <-
        protein_data_gene_upstream[i, ]$End
      message(
        paste0(
          "The protein domain ",
          protein_data_gene_upstream[i, ]$Domain_name_abbreviation,
          " from the upstream gene (",
          fusion@gene_upstream@name,
          ") is retained in the fusion transcript!"
        )
      )
    } else if (
      protein_data_gene_upstream[i, ]$Start >=
        min(start(exons_part_of_fusion_downshifted_a)) &
      protein_data_gene_upstream[i, ]$Start <
        max(end(exons_part_of_fusion_downshifted_a))
    ) {
      # Is it partly inside?
      protein_data_gene_upstream[i, ]$protein_domain_location <-
        "insideButBroken"
      protein_data_gene_upstream[i, ]$plot_start <-
        protein_data_gene_upstream[i, ]$Start
      protein_data_gene_upstream[i, ]$plot_end <-
        fusion_breakpoint_in_downshifted_transcript_a
    } else {
      protein_data_gene_upstream[i, ]$protein_domain_location <- "outside"
    }
  }
  if (
    nrow(
      dplyr::filter(
        protein_data_gene_upstream,
        protein_domain_location == "inside"
      )
    ) == 0
  ) {
    message("No protein domains are retained in the upstream gene.")
  }
  # b)
  if (fusion@inframe ||
      plot_downstream_protein_domains_if_fusion_is_out_of_frame
    ) {
    # We only consider downstream protein domains if the reading frame is kept
    # in the fusion, or if the user specifically asks for it.
    if (!fusion@inframe &&
        plot_downstream_protein_domains_if_fusion_is_out_of_frame
      ) {
      message(
        paste0(
          "No protein domains are retained in the downstream gene, because ",
          "the fusion is not in frame. But since ",
          "plot_downstream_protein_domains_if_fusion_is_out_of_frame = TRUE,",
          " we'll plot the protein domains in the downstream gene anyway.."
        )
      )
    }
    for (i in 1:nrow(protein_data_gene_downstream)) {
      if (
        protein_data_gene_downstream[i, ]$Start >=
          min(start(exons_part_of_fusion_downshifted_b)) &
        protein_data_gene_downstream[i, ]$End <=
          max(end(exons_part_of_fusion_downshifted_b))
      ) {
        # Is it "inside"?
        protein_data_gene_downstream[i, ]$protein_domain_location <- "inside"
        # When we connect the exons from geneA that are part of the fusion with
        # the exons from geneB that are part of the fusion, then the coordinates
        # we'll plot the protein domain on will shift to the right (as many
        # steps as the total length of the exons that are included from geneA).
        # That is why we add sum(width(exons_part_of_fusion_downshifted_a)) -
        # the total length of the exons that are included from geneA - below:
        protein_data_gene_downstream[i, ]$plot_start <-
          protein_data_gene_downstream[i, ]$Start +
            sum(width(exons_part_of_fusion_downshifted_a))
        protein_data_gene_downstream[i, ]$plot_end <-
          protein_data_gene_downstream[i, ]$End +
            sum(width(exons_part_of_fusion_downshifted_a))
        message(
          paste0(
            "The protein domain ",
            protein_data_gene_downstream[i, ]$Domain_name_abbreviation,
            " from the downstream gene (",
            fusion@gene_downstream@name,
            ") is retained in the fusion transcript!"
          )
        )
      } else if (
        protein_data_gene_downstream[i, ]$Start >=
          min(start(exons_part_of_fusion_downshifted_b)) &
        protein_data_gene_downstream[i, ]$Start <
          max(end(exons_part_of_fusion_downshifted_b))
      ) {
        # Is it partly inside?
        protein_data_gene_downstream[i, ]$protein_domain_location <-
          "insideButBroken"
        protein_data_gene_downstream[i, ]$plot_start <-
          protein_data_gene_downstream[i, ]$Start +
            sum(width(exons_part_of_fusion_downshifted_a))
        protein_data_gene_downstream[i, ]$plot_end <-
          fusion_breakpoint_in_downshifted_transcript_b
      } else {
        protein_data_gene_downstream[i, ]$protein_domain_location <- "outside"
      }
    }
    if (
      nrow(
        dplyr::filter(
          protein_data_gene_downstream,
          protein_domain_location == "inside"
        )
      ) == 0
    ) {
      message("No protein domains are retained in the downstream gene.")
    }
  } else {
    message(paste0("The fusion is out of frame, so no protein domains are ",
                   "retained in the downstream gene."))
  }

  # Reverse order of exons if on minus strand
  if (fusion@gene_upstream@strand == "-") {
    fusion_exons_upstream <- rev(exons_part_of_fusion_a)
  } else {
    fusion_exons_upstream <- exons_part_of_fusion_a
  }
  if (fusion@gene_downstream@strand == "-") {
    fusion_exons_downstream <- rev(exons_part_of_fusion_b)
  } else {
    fusion_exons_downstream <- exons_part_of_fusion_b
  }

  # Build fusion transcript
  fusion_transcript <- append(
    fusion_exons_upstream,
    fusion_exons_downstream
  )

  # If we've got a bamfile, calculate coverage
  if (!is.null(bamfile)) {
    # Get coverage
    cov <- coverage(
      bamfile,
      param = Rsamtools::ScanBamParam(
        which = fusion_transcript
      )
    )    ####  This must be unshifted!
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
  }

  # Downshift the fusion transcript (remove introns)
  fusion_transcript <- down_shift(fusion_transcript)

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
    showTitle = FALSE, # hide title
    thinBoxFeature = c(
      "5utr_excludedA",
      "5utr_includedA",
      "3utr_excludedA",
      "3utr_includedA",

      "5utr_excludedB",
      "5utr_includedB",
      "3utr_excludedB",
      "3utr_includedB"
    )
  )

  # Color exons
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
  Gviz::feature(tr_track)[
    Gviz::symbol(tr_track) == fusion@gene_upstream@name
  ] <- "protein_coding_includedA"
  Gviz::feature(tr_track)[
    Gviz::symbol(tr_track) == fusion@gene_upstream@name &
      mcols(fusion_transcript)$feature == "5utr"
  ] <- "5utr_includedA"
  Gviz::feature(tr_track)[
    Gviz::symbol(tr_track) == fusion@gene_upstream@name &
      mcols(fusion_transcript)$feature == "3utr"
  ] <- "3utr_includedA"
  Gviz::feature(tr_track)[
    Gviz::symbol(tr_track) == fusion@gene_downstream@name
  ] <- "protein_coding_includedB"
  Gviz::feature(tr_track)[
    Gviz::symbol(tr_track) == fusion@gene_downstream@name &
      mcols(fusion_transcript)$feature == "5utr"
  ] <- "5utr_includedB"
  Gviz::feature(tr_track)[
    Gviz::symbol(tr_track) == fusion@gene_downstream@name &
      mcols(fusion_transcript)$feature == "3utr"
  ] <- "3utr_includedB"
  Gviz::displayPars(tr_track) <- interestcolor

  # Create axis track
  axis_track <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axis_track) <- list(
    col = "black", # line color
    fontcolor = "black",
    fontsize =  25,
    lwd = 1.5,
    showId = TRUE
  )

  # What are the start and end positions of the protein domain annotations?
  protein_domain_annotations_start <- c(
    dplyr::filter(
      protein_data_gene_upstream,
      protein_domain_location == "inside"
    )$plot_start,
    dplyr::filter(
      protein_data_gene_downstream,
      protein_domain_location == "inside"
    )$plot_start
  )
  protein_domain_annotations_end <- c(
    dplyr::filter(
      protein_data_gene_upstream,
      protein_domain_location == "inside"
    )$plot_start,
    dplyr::filter(
      protein_data_gene_downstream,
      protein_domain_location == "inside"
    )$plot_end
  )
  protein_domain_annotations_names <- c(
    dplyr::filter(
      protein_data_gene_upstream,
      protein_domain_location == "inside"
    )$plot_start,
    dplyr::filter(
      protein_data_gene_downstream,
      protein_domain_location == "inside"
    )$Domain_name_abbreviation
  )
  # Create protein domain track
  annotation_track <- Gviz::AnnotationTrack(
    start = protein_domain_annotations_start,
    end = protein_domain_annotations_end,
    chromosome = "chrNA",
    id = protein_domain_annotations_names,
    genome = "hg19")
  # Set display parameters
  Gviz::displayPars(annotation_track) <- list(
    background.panel = "transparent", # background color of the content panel
    background.title = "transparent", # background color for the title panels
    col.title = "black",
    showTitle = FALSE, # show title
    showId = FALSE, # show ids
    col = "lightgray", # border color
    showTitle = FALSE, # hide left title
    featureAnnotation = "id",
    fontcolor.feature = "black"
  )

  # If we've got a bamfile, plot with coverage
  if (!is.null(bamfile)) {
    Gviz::plotTracks(
      list(d_track, tr_track, annotation_track, axis_track),
      main = paste(
        fusion@gene_upstream@name,
        fusion@gene_downstream@name,
        sep = ":"
      ),
      chromosome = "chrNA")
  } else {
    Gviz::plotTracks(
      list(tr_track, annotation_track, axis_track),
      main = paste(
        fusion@gene_upstream@name,
        fusion@gene_downstream@name,
        sep = ":"
      ),
      chromosome = "chrNA")
  }

}

.validate_plot_fusion_transcript_with_protein_domain_params <- function(
  fusion,
  edb,
  bamfile,
  bedfile,
  gene_upstream_transcript,
  gene_downstream_transcript,
  plot_downstream_protein_domains_if_fusion_is_out_of_frame
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()

  # Check parameters
  argument_checker <- .is_fusion_valid(argument_checker, fusion)
  argument_checker <- .is_edb_valid(argument_checker, edb, fusion)
  argument_checker <- .is_bamfile_valid(argument_checker, bamfile)
  argument_checker <- .is_bedfile_valid(argument_checker, bedfile)
  argument_checker <- .is_character_parameter_valid(
    argument_checker,
    gene_upstream_transcript,
    "gene_upstream_transcript")
  argument_checker <- .is_character_parameter_valid(
    argument_checker,
    gene_downstream_transcript,
    "gene_downstream_transcript")
  argument_checker <- .is_parameter_boolean(
    argument_checker,
    plot_downstream_protein_domains_if_fusion_is_out_of_frame,
    "plot_downstream_protein_domains_if_fusion_is_out_of_frame")

  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)
}
