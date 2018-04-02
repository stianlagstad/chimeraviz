#' Import results from a SOAPfuse run into a list of Fusion objects.
#'
#' A function that imports the results from a SOAPfuse run, typically from
#' a final.Fusion.specific.for.genes file, into a list of Fusion objects.
#'
#' @param filename Filename for the SOAPfuse
#' final-list-candidate-fusion-genes.txt results file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' soapfuse833ke <- system.file(
#'   "extdata",
#'   "soapfuse_833ke_final.Fusion.specific.for.genes",
#'   package = "chimeraviz")
#' fusions <- import_soapfuse(soapfuse833ke, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
import_soapfuse <- function (filename, genome_version, limit) {

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
      col_types_soapfuse <- readr::cols_only(
        "up_gene" = readr::col_character(),
        "up_chr" = readr::col_character(),
        "up_strand" = readr::col_character(),
        "up_Genome_pos" = readr::col_integer(),
        "dw_gene" = readr::col_character(),
        "dw_chr" = readr::col_character(),
        "dw_strand" = readr::col_character(),
        "dw_Genome_pos" = readr::col_integer(),
        "Span_reads_num" = readr::col_integer(),
        "Junc_reads_num" = readr::col_integer(),
        "down_fusion_part_frame-shift_or_not" = readr::col_character()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_soapfuse
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_soapfuse,
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
  fusion_tool          <- "soapfuse"
  spanning_reads_count <- NA
  split_reads_count    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import soapfuse-specific fields
    fusion_tool_specific_data <- list()

    # Set id for this fusion. Soapfuse does not have any unique identifier
    # for each row of the results file, so just use our counter.
    id <- as.character(i)

    # Is the downstream fusion partner in-frame?
    if (is.na(report[[i, "down_fusion_part_frame-shift_or_not"]])) {
      inframe <- NA
    } else if (
      report[[i, "down_fusion_part_frame-shift_or_not"]] == "inframe-shift"
    ) {
      inframe <- TRUE
    } else {
      inframe <- FALSE
    }

    # Chromosome names
    chromosome_upstream <- report[[i, "up_chr"]]
    chromosome_downstream <- report[[i, "dw_chr"]]

    # Breakpoints
    breakpoint_upstream <- report[[i, "up_Genome_pos"]]
    breakpoint_downstream <- report[[i, "dw_Genome_pos"]]

    # Strand
    strand_upstream <- report[[i, "up_strand"]]
    strand_downstream <- report[[i, "dw_strand"]]

    # Number of supporting reads
    split_reads_count <- report[[i, "Span_reads_num"]]
    spanning_reads_count <- report[[i, "Junc_reads_num"]]

    # SOAPfuse doesn't provide the fusion sequence
    junction_sequence_upstream <- Biostrings::DNAString()
    junction_sequence_downstream <- Biostrings::DNAString()

    # Gene names
    name_upstream <- report[[i, "up_gene"]]
    name_downstream <- report[[i, "dw_gene"]]

    # Ensembl ids
    ensembl_id_upstream <- NA_character_
    ensembl_id_downstream <- NA_character_

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
