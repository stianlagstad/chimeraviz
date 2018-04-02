#' Import results from a EricScript run into a list of Fusion objects.
#'
#' A function that imports the results from a EricScript run into a list of
#' Fusion objects.
#'
#' @param filename Filename for the EricScript results file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' ericscriptData <- system.file(
#'   "extdata",
#'   "ericscript_SRR1657556.results.total.tsv",
#'   package = "chimeraviz")
#' fusions <- import_ericscript(ericscriptData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
import_ericscript <- function (filename, genome_version, limit) {

  # Is the genome version valid?
  valid_genomes <- c("hg19", "hg38", "mm10")
  if (is.na(match(tolower(genome_version), tolower(valid_genomes)))) {
    stop("Invalid genome version given")
  }

  # If the limit is set, is the value valid?
  if (missing(limit) == FALSE) {
    if (is.numeric(limit) == FALSE || limit <= 0) {
      stop("limit must be a numeric value bigger than 0")
    }
  }

  # Try to read the fusion report
  report <- withCallingHandlers({
      col_types_ericscript <- readr::cols_only(
        "GeneName1" = col_character(),
        "GeneName2" = col_character(),
        "chr1" = col_character(),
        "Breakpoint1" = col_character(),
        "strand1" = col_character(),
        "chr2" = col_character(),
        "Breakpoint2" = col_character(),
        "strand2" = col_character(),
        "EnsemblGene1" = col_character(),
        "EnsemblGene2" = col_character(),
        "crossingreads" = col_integer(),
        "spanningreads" = col_integer(),
        "mean.insertsize" = col_skip(),
        "homology" = col_skip(),
        "fusiontype" = col_skip(),
        "Blacklist" = col_skip(),
        "InfoGene1" = col_skip(),
        "InfoGene2" = col_skip(),
        "JunctionSequence" = col_character(),
        "GeneExpr1" = col_skip(),
        "GeneExpr2" = col_skip(),
        "GeneExpr_Fused" = col_skip(),
        "ES" = col_skip(),
        "GJS" = col_skip(),
        "US" = col_skip(),
        "EricScore" = col_number()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_ericscript
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_ericscript,
          n_max = limit
        )
      }
    },
    error = function(cond) {
      message(paste0("Reading ", filename, " caused an error: ", cond[[1]]))
      stop(cond)
    },
    warning = function(cond) {
      # Begin Exclude Linting
      #message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      # End Exclude Linting
    }
  )

  # Set variables
  id                   <- NA
  inframe              <- NA
  fusion_tool          <- "ericscript"
  spanning_reads_count <- NA
  split_reads_count    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import ericscript-specific fields
    fusion_tool_specific_data <- list()
    fusion_tool_specific_data[["EricScore"]] <- report[[i, "EricScore"]]

    # Set id for this fusion. Soapfuse does not have any unique identifier
    # for each row of the results file, so just use our counter.
    id <- as.character(i)

    # Chromosome names
    chromosome_upstream <- paste("chr", report[[i, "chr1"]], sep = "")
    chromosome_downstream <- paste("chr", report[[i, "chr2"]], sep = "")

    # Breakpoints
    # EricScript doesn't always report fusion breakpoints. If missing, then set
    # breakpoint to -1
    breakpoint_upstream <- tryCatch(
      as.numeric(report[[i, "Breakpoint1"]]),
      warning = function(w) {
        -1
      }
    )
    breakpoint_downstream <- tryCatch(
      as.numeric(report[[i, "Breakpoint2"]]),
      warning = function(w) {
        -1
      }
    )

    # Strand
    strand_upstream <- report[[i, "strand1"]]
    strand_downstream <- report[[i, "strand2"]]

    # Number of supporting reads
    split_reads_count <- report[[i, "crossingreads"]]
    spanning_reads_count <- report[[i, "spanningreads"]]

    # Get the fusion sequence
    # EricScript outputs this as sequencebeforejunctionSEQUENCEAFTERJUNCTION, so
    # first locate the first upper cased letter:
    split_at <-
      regexpr(pattern = "([[:upper:]])", report[[i, "JunctionSequence"]])[1]
    # Then we can extract each piece:
    junction_sequence_upstream <- Biostrings::DNAString(
      substr(
        report[[i, "JunctionSequence"]],
        1,
        split_at - 1)
    )
    junction_sequence_downstream <- Biostrings::DNAString(
      substr(
        report[[i, "JunctionSequence"]],
        split_at,
        nchar(report[[i, "JunctionSequence"]]))
    )

    # Gene names
    name_upstream <- report[[i, "GeneName1"]]
    name_downstream <- report[[i, "GeneName2"]]

    # Ensembl ids
    ensembl_id_upstream <- report[[i, "EnsemblGene1"]]
    ensembl_id_downstream <- report[[i, "EnsemblGene2"]]

    # PartnerGene objects
    gene_upstream <- new(
      Class = "PartnerGene",
      name = name_upstream,
      ensembl_id = ensembl_id_upstream,
      chromosome = chromosome_upstream,
      breakpoint = breakpoint_upstream,
      strand = strand_upstream,
      junction_sequence = junction_sequence_upstream,
      transcripts = GenomicRanges::GRangesList()
    )
    gene_downstream <- new(
      Class = "PartnerGene",
      name = name_downstream,
      ensembl_id = ensembl_id_downstream,
      chromosome = chromosome_downstream,
      breakpoint = breakpoint_downstream,
      strand = strand_downstream,
      junction_sequence = junction_sequence_downstream,
      transcripts = GenomicRanges::GRangesList()
    )

    fusion_list[[i]] <- new(
      Class = "Fusion",
      id = id,
      fusion_tool = fusion_tool,
      genome_version = genome_version,
      spanning_reads_count = spanning_reads_count,
      split_reads_count = split_reads_count,
      fusion_reads_alignment = Gviz::AlignmentsTrack(),
      gene_upstream = gene_upstream,
      gene_downstream = gene_downstream,
      inframe = inframe,
      fusion_tool_specific_data = fusion_tool_specific_data
    )
  }

  # Return the list of Fusion objects
  fusion_list
}
