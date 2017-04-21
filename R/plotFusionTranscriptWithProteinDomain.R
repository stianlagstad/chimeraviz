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
#' @param geneAtranscript The transcript id for the upstream gene.
#' @param geneBtranscript The transcript id for the downstream gene.
#' @param plotDownstreamProteinDomainsIfFusionIsOutOfFrame Setting this to true
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
#' fusions <- importDefuse(defuse833ke, "hg19", 1)
#' fusion <- getFusionById(fusions, 5267)
#' # Select transcripts
#' geneAtranscript <- "ENST00000434290"
#' geneBtranscript <- "ENST00000370031"
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
#' plotFusionTranscriptWithProteinDomain(
#'   fusion = fusion,
#'   edb = edb,
#'   bamfile = bamfile5267,
#'   bedfile = bedfile,
#'   geneAtranscript = geneAtranscript,
#'   geneBtranscript = geneBtranscript,
#'   plotDownstreamProteinDomainsIfFusionIsOutOfFrame = TRUE)
#' # Close device
#' dev.off()
#'
#' @export
#'
plotFusionTranscriptWithProteinDomain <- function(
  fusion,
  edb = NULL,
  bamfile = NULL,
  bedfile = NULL,
  geneAtranscript = "",
  geneBtranscript = "",
  plotDownstreamProteinDomainsIfFusionIsOutOfFrame = FALSE) {

  # Validate that we're given a bedfile that exists
  if (class(bedfile) != "character") {
    stop("bedfile argument must be a character")
  }
  if (file.exists(bedfile) == FALSE) {
    stop("Invalid bedfile: File doesn't exist.")
  }
  # Validate geneAtranscript and geneBtranscript
  if (class(geneAtranscript) != "character") {
    stop("geneAtranscript argument must be a character")
  }
  if (class(geneBtranscript) != "character") {
    stop("geneBtranscript argument must be a character")
  }
  if (geneAtranscript == "") {
    stop("geneAtranscript argument must not be blank")
  }
  if (geneBtranscript == "") {
    stop("geneBtranscript argument must not be blank")
  }
  # Validate parameters and fetch transcripts from database
  fusion <- .validateFusionPlotParams(
    fusion,
    edb,
    bamfile)
  # Validate that the transcript ids are found in the edb
  if (length(fusion@geneA@transcripts[names(fusion@geneA@transcripts) == geneAtranscript]) == 0) {
    stop(paste0(geneAtranscript, " was not found in among the transcripts for the upstream gene partner."))
  }
  if (length(fusion@geneB@transcripts[names(fusion@geneB@transcripts) == geneBtranscript]) == 0) {
    stop(paste0(geneBtranscript, " was not found in among the transcripts for the downstream gene partner."))
  }
  # Extract the transcripts
  transcriptA <- unlist(fusion@geneA@transcripts[names(fusion@geneA@transcripts) == geneAtranscript])
  transcriptB <- unlist(fusion@geneB@transcripts[names(fusion@geneB@transcripts) == geneBtranscript])
  # Check that there's at least one transcript that has the fusion breakpoint
  # within the transcript.
  .checkThatBreakpointsAreWithinTranscripts(fusion, transcriptA, transcriptB)
  # Check that the transcripts have a breakpoint exon
  .checkThatTranscriptsHaveBreakpointExons(fusion, transcriptA, transcriptB)
  # Add transcript metadata column to make Gviz group correctly
  mcols(transcriptA)$transcript <- fusion@geneA@name
  mcols(transcriptB)$transcript <- fusion@geneB@name
  # Add symbol metadata column to make Gviz show ids
  mcols(transcriptA)$symbol <- fusion@geneA@name
  mcols(transcriptB)$symbol <- fusion@geneB@name

  # Downshift transcripts (remove introns from the GRanges)
  transcriptAdownshifted <- downShift(transcriptA)
  transcriptBdownshifted <- downShift(transcriptB)

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
  #    - Exons part of the fusion transcript are those with start/end positions before breakpoint
  #    - Fusion breakpoint in the downshifted transcript is at the end
  #      coordinate of exon_id_at_fusion_breakpoint_b
  #
  if (fusion@geneA@strand == "+") {

    exon_ids_part_of_fusion_from_transcript_a <- mcols(transcriptA[start(ranges(transcriptA)) < fusion@geneA@breakpoint & end(ranges(transcriptA)) <= fusion@geneA@breakpoint])$exon_id
    exons_part_of_fusion_a <- transcriptA[start(ranges(transcriptA)) < fusion@geneA@breakpoint & end(ranges(transcriptA)) <= fusion@geneA@breakpoint]
    exon_id_at_fusion_breakpoint_a <- mcols(transcriptA[end(ranges(transcriptA)) == fusion@geneA@breakpoint])$exon_id

    fusion_breakpoint_in_downshifted_transcript_a <- max(end(transcriptAdownshifted[mcols(transcriptAdownshifted)$exon_id == exon_id_at_fusion_breakpoint_a]))
    exons_part_of_fusion_downshifted_a <- transcriptAdownshifted[start(ranges(transcriptAdownshifted)) < fusion_breakpoint_in_downshifted_transcript_a & end(ranges(transcriptAdownshifted)) <= fusion_breakpoint_in_downshifted_transcript_a]

  } else {

    exon_ids_part_of_fusion_from_transcript_a <- mcols(transcriptA[start(ranges(transcriptA)) >= fusion@geneA@breakpoint & end(ranges(transcriptA)) > fusion@geneA@breakpoint])$exon_id
    exons_part_of_fusion_a <- transcriptA[start(ranges(transcriptA)) >= fusion@geneA@breakpoint & end(ranges(transcriptA)) > fusion@geneA@breakpoint]
    exon_id_at_fusion_breakpoint_a <- mcols(transcriptA[start(ranges(transcriptA)) == fusion@geneA@breakpoint])$exon_id

    fusion_breakpoint_in_downshifted_transcript_a <- min(start(transcriptAdownshifted[mcols(transcriptAdownshifted)$exon_id == exon_id_at_fusion_breakpoint_a]))
    exons_part_of_fusion_downshifted_a <- transcriptAdownshifted[start(ranges(transcriptAdownshifted)) >= fusion_breakpoint_in_downshifted_transcript_a & end(ranges(transcriptAdownshifted)) > fusion_breakpoint_in_downshifted_transcript_a]

  }
  if (fusion@geneB@strand == "+") {

    exon_ids_part_of_fusion_from_transcript_b <- mcols(transcriptB[start(ranges(transcriptB)) >= fusion@geneB@breakpoint & end(ranges(transcriptB)) > fusion@geneB@breakpoint])$exon_id
    exons_part_of_fusion_b <- transcriptB[start(ranges(transcriptB)) >= fusion@geneB@breakpoint & end(ranges(transcriptB)) > fusion@geneB@breakpoint]
    exon_id_at_fusion_breakpoint_b <- mcols(transcriptB[start(ranges(transcriptB)) == fusion@geneB@breakpoint])$exon_id

    fusion_breakpoint_in_downshifted_transcript_b <- min(start(transcriptBdownshifted[mcols(transcriptBdownshifted)$exon_id == exon_id_at_fusion_breakpoint_b]))
    exons_part_of_fusion_downshifted_b <- transcriptBdownshifted[start(ranges(transcriptBdownshifted)) >= fusion_breakpoint_in_downshifted_transcript_b & end(ranges(transcriptBdownshifted)) > fusion_breakpoint_in_downshifted_transcript_b]

  } else {

    exon_ids_part_of_fusion_from_transcript_b <- mcols(transcriptB[start(ranges(transcriptB)) < fusion@geneB@breakpoint & end(ranges(transcriptB)) <= fusion@geneB@breakpoint])$exon_id
    exons_part_of_fusion_b <- transcriptB[start(ranges(transcriptB)) < fusion@geneB@breakpoint & end(ranges(transcriptB)) <= fusion@geneB@breakpoint]
    exon_id_at_fusion_breakpoint_b <- mcols(transcriptB[end(ranges(transcriptB)) == fusion@geneB@breakpoint])$exon_id

    fusion_breakpoint_in_downshifted_transcript_b <- max(end(transcriptBdownshifted[mcols(transcriptBdownshifted)$exon_id == exon_id_at_fusion_breakpoint_b]))
    exons_part_of_fusion_downshifted_b <- transcriptBdownshifted[start(ranges(transcriptBdownshifted)) < fusion_breakpoint_in_downshifted_transcript_b & end(ranges(transcriptBdownshifted)) <= fusion_breakpoint_in_downshifted_transcript_b]

  }

  # Load protein domain data
  # Specify which columns from the bedfile we're going to read
  col_types_protein_data <- readr::cols_only(
    "Transcript_id" = readr::col_character(),
    "Pfam_id" = readr::col_character(),
    "Start" = readr::col_integer(),
    "End" = readr::col_integer(),
    "Domain_name_abbreviation" = readr::col_character(),
    "Domain_name_full" = readr::col_character()
  )
  # Read the bedfile
  proteinData <- readr::read_tsv(
    file = bedfile,
    col_types = col_types_protein_data
  )
  # Validate that the transcript ids are found in the bedfile
  if (nrow(proteinData[proteinData$Transcript_id == geneAtranscript, ]) == 0) {
    stop(paste0("Found no transcripts matching ", geneAtranscript, " in the edb database."))
  }
  if (nrow(proteinData[proteinData$Transcript_id == geneBtranscript, ]) == 0) {
    stop(paste0("Found no transcripts matching ", geneBtranscript, " in the edb database."))
  }
  # Extract the rows we need
  proteinDataGeneA <- proteinData[proteinData$Transcript_id == geneAtranscript, ]
  proteinDataGeneB <- proteinData[proteinData$Transcript_id == geneBtranscript, ]

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
  # 1: geneA is on the + strand
  #    - The protein domain is "inside" if: protein_start_coordinate >= min(start(exons_part_of_fusion_downshifted_a)) & protein_end_coordinate <= max(end(exons_part_of_fusion_downshifted_a))
  #    - The protein domain is "insideButBroken" if: protein_start_coordinate >= min(start(exons_part_of_fusion_downshifted_a)) & protein_start_coordinate < max(end(exons_part_of_fusion_downshifted_a))
  #
  # 2: geneA is on the - strand
  #    - The protein domain is "inside" if: protein_start_coordinate <= max(end(exons_part_of_fusion_downshifted_a)) & protein_end_coordinate >= min(start(exons_part_of_fusion_downshifted_a))
  #    - The protein domain is "insideButBroken" if: protein_start_coordinate <= max(end(exons_part_of_fusion_downshifted_a)) & protein_start_coordinate > min(start(exons_part_of_fusion_downshifted_a))
  #
  # 3: geneB is on the + strand
  #    - The protein domain is "inside" if: protein_start_coordinate >= min(start(exons_part_of_fusion_downshifted_b)) & protein_end_coordinate <= max(end(exons_part_of_fusion_downshifted_b))
  #    - The protein domain is "insideButBroken" if: protein_start_coordinate >= min(start(exons_part_of_fusion_downshifted_b)) & protein_start_coordinate < max(end(exons_part_of_fusion_downshifted_b))
  #
  # 4: geneB is on the - strand
  #    - The protein domain is "inside" if: protein_start_coordinate <= max(end(exons_part_of_fusion_downshifted_b)) & protein_end_coordinate >= min(start(exons_part_of_fusion_downshifted_b))
  #    - The protein domain is "insideButBroken" if: protein_start_coordinate <= max(end(exons_part_of_fusion_downshifted_b)) & protein_start_coordinate > min(start(exons_part_of_fusion_downshifted_b))

  # First, if the gene is on the minus strand, then we have to "shift" the Start
  # end End coordinates we have for the protein domain:
  # if (fusion@geneA@strand == "-") {
  #   new_start <- sum(width(transcriptAdownshifted)) - proteinDataGeneA$End
  #   new_end <- sum(width(transcriptAdownshifted)) - proteinDataGeneA$Start
  #   proteinDataGeneA$Start <- new_start
  #   proteinDataGeneA$End <- new_end
  # }
  # if (fusion@geneB@strand == "-") {
  #   new_start <- sum(width(transcriptBdownshifted)) - proteinDataGeneB$End
  #   new_end <- sum(width(transcriptBdownshifted)) - proteinDataGeneB$Start
  #   proteinDataGeneB$Start <- new_start
  #   proteinDataGeneB$End <- new_end
  # }

  # Go through all protein domains to classify them.
  # First extend the proteinDataGeneA and proteinDataGeneB data frames:
  proteinDataGeneA$protein_domain_location <- "none"
  proteinDataGeneA$plot_start <- "-"
  proteinDataGeneA$plot_end <- "-"
  proteinDataGeneB$protein_domain_location <- "none"
  proteinDataGeneB$plot_start <- "-"
  proteinDataGeneB$plot_end <- "-"
  # a)
  for (i in 1:nrow(proteinDataGeneA)) {
    if (proteinDataGeneA[i,]$Start >= min(start(exons_part_of_fusion_downshifted_a)) & proteinDataGeneA[i,]$End <= max(end(exons_part_of_fusion_downshifted_a))) {
      # Is it "inside"?
      proteinDataGeneA[i,]$protein_domain_location <- "inside"
      proteinDataGeneA[i,]$plot_start <- proteinDataGeneA[i,]$Start
      proteinDataGeneA[i,]$plot_end <- proteinDataGeneA[i,]$End
      message(paste0("The protein domain ", proteinDataGeneA[i,]$Domain_name_abbreviation, " from the upstream gene (", fusion@geneA@name, ") is retained in the fusion transcript!"))
    } else if (proteinDataGeneA[i,]$Start >= min(start(exons_part_of_fusion_downshifted_a)) & proteinDataGeneA[i,]$Start < max(end(exons_part_of_fusion_downshifted_a))) {
      # Is it partly inside?
      proteinDataGeneA[i,]$protein_domain_location <- "insideButBroken"
      proteinDataGeneA[i,]$plot_start <- proteinDataGeneA[i,]$Start
      proteinDataGeneA[i,]$plot_end <- fusion_breakpoint_in_downshifted_transcript_a
      # message(paste0("The protein domain ", proteinDataGeneA[i,]$Domain_name_abbreviation, " from the upstream gene (", fusion@geneA@name, ") is partly retained in the fusion transcript."))
    } else {
      proteinDataGeneA[i,]$protein_domain_location <- "outside"
    }
  }
  if (nrow(dplyr::filter(proteinDataGeneA, protein_domain_location == "inside")) == 0) {
    message("No protein domains are retained in the upstream gene.")
  }
  # b)
  if (fusion@inframe || plotDownstreamProteinDomainsIfFusionIsOutOfFrame) {
    # We only consider downstream protein domains if the reading frame is kept
    # in the fusion, or if the user specifically asks for it.
    if (!fusion@inframe && plotDownstreamProteinDomainsIfFusionIsOutOfFrame) {
      message(paste0("The fusion is out of frame, but since ",
                     "plotDownstreamProteinDomainsIfFusionIsOutOfFrame = TRUE,",
                     " we'll consider them anyway.."))
    }
    for (i in 1:nrow(proteinDataGeneB)) {
      if (proteinDataGeneB[i,]$Start >= min(start(exons_part_of_fusion_downshifted_b)) & proteinDataGeneB[i,]$End <= max(end(exons_part_of_fusion_downshifted_b))) {
        # Is it "inside"?
        proteinDataGeneB[i,]$protein_domain_location <- "inside"
        # When we connect the exons from geneA that are part of the fusion with
        # the exons from geneB that are part of the fusion, then the coordinates
        # we'll plot the protein domain on will shift to the right (as many
        # steps as the total length of the exons that are included from geneA).
        # That is why we add sum(width(exons_part_of_fusion_downshifted_a)) -
        # the total length of the exons that are included from geneA - below:
        proteinDataGeneB[i,]$plot_start <- proteinDataGeneB[i,]$Start + sum(width(exons_part_of_fusion_downshifted_a))
        proteinDataGeneB[i,]$plot_end <- proteinDataGeneB[i,]$End + sum(width(exons_part_of_fusion_downshifted_a))
        message(paste0("The protein domain ", proteinDataGeneB[i,]$Domain_name_abbreviation, " from the downstream gene (", fusion@geneB@name, ") is retained in the fusion transcript!"))
      } else if (proteinDataGeneB[i,]$Start >= min(start(exons_part_of_fusion_downshifted_b)) & proteinDataGeneB[i,]$Start < max(end(exons_part_of_fusion_downshifted_b))) {
        # Is it partly inside?
        proteinDataGeneB[i,]$protein_domain_location <- "insideButBroken"
        proteinDataGeneB[i,]$plot_start <- proteinDataGeneB[i,]$Start + sum(width(exons_part_of_fusion_downshifted_a))
        proteinDataGeneB[i,]$plot_end <- fusion_breakpoint_in_downshifted_transcript_b
        # message(paste0("The protein domain ", proteinDataGeneB[i,]$Domain_name_abbreviation, " from the downstream gene (", fusion@geneB@name, ") is partly retained in the fusion transcript."))
      } else {
        proteinDataGeneB[i,]$protein_domain_location <- "outside"
      }
    }
    if (nrow(dplyr::filter(proteinDataGeneB, protein_domain_location == "inside")) == 0) {
      message("No protein domains are retained in the downstream gene.")
    }
  } else {
    message(paste0("The fusion is out of frame, so no protein domains are ",
                   "retained in the downstream gene."))
  }

  # Reverse order of exons if on minus strand
  if (fusion@geneA@strand == "-") {
    fusionExonsA <- rev(exons_part_of_fusion_a)
  } else {
    fusionExonsA <- exons_part_of_fusion_a
  }
  if (fusion@geneB@strand == "-") {
    fusionExonsB <- rev(exons_part_of_fusion_b)
  } else {
    fusionExonsB <- exons_part_of_fusion_b
  }

  # Build fusion transcript
  fusionTranscript <- append(
    fusionExonsA,
    fusionExonsB
  )

  # If we've got a bamfile, calculate coverage
  if (!is.null(bamfile)) {
    # Get coverage
    cov <- coverage(bamfile, param = Rsamtools::ScanBamParam(which = fusionTranscript))    ####  This must be unshifted!
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

  # Downshift the fusion transcript (remove introns)
  fusionTranscript <- downShift(fusionTranscript)

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
  Gviz::feature(trTrack)[Gviz::symbol(trTrack) == fusion@geneA@name] <- "protein_coding_includedA"
  Gviz::feature(trTrack)[Gviz::symbol(trTrack) == fusion@geneA@name & mcols(fusionTranscript)$feature == "5utr"] <- "5utr_includedA"
  Gviz::feature(trTrack)[Gviz::symbol(trTrack) == fusion@geneA@name & mcols(fusionTranscript)$feature == "3utr"] <- "3utr_includedA"
  Gviz::feature(trTrack)[Gviz::symbol(trTrack) == fusion@geneB@name] <- "protein_coding_includedB"
  Gviz::feature(trTrack)[Gviz::symbol(trTrack) == fusion@geneB@name & mcols(fusionTranscript)$feature == "5utr"] <- "5utr_includedB"
  Gviz::feature(trTrack)[Gviz::symbol(trTrack) == fusion@geneB@name & mcols(fusionTranscript)$feature == "3utr"] <- "3utr_includedB"
  Gviz::displayPars(trTrack) <- interestcolor

  # Create axis track
  axisTrack <- Gviz::GenomeAxisTrack()
  Gviz::displayPars(axisTrack) <- list(
    col = "black", # line color
    fontcolor = "black",
    fontsize= 25,
    lwd = 1.5,
    showId = TRUE
  )

  # What are the start and end positions of the protein domain annotations?
  protein_domain_annotations_start <- c(dplyr::filter(proteinDataGeneA, protein_domain_location == "inside")$plot_start, dplyr::filter(proteinDataGeneB, protein_domain_location == "inside")$plot_start)
  protein_domain_annotations_end <- c(dplyr::filter(proteinDataGeneA, protein_domain_location == "inside")$plot_start, dplyr::filter(proteinDataGeneB, protein_domain_location == "inside")$plot_end)
  protein_domain_annotations_names <- c(dplyr::filter(proteinDataGeneA, protein_domain_location == "inside")$plot_start, dplyr::filter(proteinDataGeneB, protein_domain_location == "inside")$Domain_name_abbreviation)
  # Create protein domain track
  annotationTrack <- Gviz::AnnotationTrack(
    start = protein_domain_annotations_start,
    end = protein_domain_annotations_end,
    chromosome = "chrNA",
    id = protein_domain_annotations_names,
    genome = "hg19")
  # Set display parameters
  Gviz::displayPars(annotationTrack) <- list(
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
      list(dTrack, trTrack, annotationTrack, axisTrack),
      main = paste(fusion@geneA@name, fusion@geneB@name, sep = ":"),
      chromosome = "chrNA")
  } else {
    Gviz::plotTracks(
      list(trTrack, annotationTrack, axisTrack),
      main = paste(fusion@geneA@name, fusion@geneB@name, sep = ":"),
      chromosome = "chrNA")
  }

}
