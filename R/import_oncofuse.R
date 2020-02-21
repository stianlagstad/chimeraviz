#' Import results from a oncofuse run into a list of Fusion objects.
#'
#' A function that imports the results from a oncofuse run, typically from a
#' results.filtered.tsv file, into a list of Fusion objects.
#' 
#' This import function was contributed by Lavinia G, ref
#' https://github.com/stianlagstad/chimeraviz/issues/47#issuecomment-409773158
#'
#' @param filename Filename for the oncofuse results .tsv file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' oncofuse833ke <- system.file(
#'   "extdata",
#'   "oncofuse.outfile",
#'   package="chimeraviz")
#' fusions <- import_oncofuse(oncofuse833ke, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
import_oncofuse <- function (filename, genome_version, limit) {

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
      col_types <- c(
        "FUSION_ID" = "character",
        "SPANNING_READS" = "integer",
        "ENCOMPASSING_READS" = "integer",
        "GENOMIC" = "character",
        "5_FPG_GENE_NAME" = "character",
        "3_FPG_GENE_NAME" = "character",
        "FPG_FRAME_DIFFERENCE" = "integer",
        "DRIVER_PROB" = "numeric"
      )
      if (missing(limit)) {
        # Read all lines
        data.table::fread(
          input = filename,
          colClasses = col_types,
          showProgress = FALSE
        )
      } else {
        # Only read up to the limit
        data.table::fread(
          input = filename,
          colClasses = col_types,
          showProgress = FALSE,
          nrows = limit
        )
      }
    },
    error = function(cond) {
      message(paste0("Reading ", filename, " caused an error: ", cond[[1]]))
      stop(cond)
    },
    warning = function(cond) {
      # Begin Exclude Linting
      message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      # End Exclude Linting
    }
  )


  # Set variables
  id                   <- NA
  inframe              <- NA
  fusion_tool          <- "oncofuse"
  spanning_reads_count <- NA
  split_reads_count    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import oncofuse-specific fields
    fusion_tool_specific_data <- list()
    fusion_tool_specific_data[["DRIVER_PROB"]] <- report[[i, "DRIVER_PROB"]]

    # Cluster id
    id <- report[[i, "FUSION_ID"]]

    # Is the downstream fusion partner in-frame?
    if (report[[i, "FPG_FRAME_DIFFERENCE"]] == "0") {
      inframe <- TRUE
    } else {
      inframe <- FALSE
    }

    # Number of supporting reads
    split_reads_count <- report[[i, "ENCOMPASSING_READS"]]
    spanning_reads_count <- report[[i, "SPANNING_READS"]]

    # Breakpoints
    genomics <- report[[i, "GENOMIC"]]
    genomic_vals <- unlist(strsplit(genomics, ">"))
    chromosome_upstream <- unlist(strsplit(genomic_vals[1], ":"))[1]
    breakpoint_upstream <- as.numeric(unlist(strsplit(genomic_vals[1], ":"))[2])
    chromosome_downstream <- unlist(strsplit(genomic_vals[2], ":"))[1]
    breakpoint_downstream <- as.numeric(
      unlist(strsplit(genomic_vals[2], ":"))[2]
    )

    # Gene names
    name_upstream <- report[[i, "5_FPG_GENE_NAME"]]
    name_downstream <- report[[i, "3_FPG_GENE_NAME"]]

    # Ensembl ids
    ensembl_id_upstream <- NA_character_
    ensembl_id_downstream <- NA_character_

    # oncofuse doesn't provide the fusion sequence
    junction_sequence_upstream <- Biostrings::DNAString()
    junction_sequence_downstream <- Biostrings::DNAString()

    # oncofuse doesn't provide the strands
    strand_upstream <- NA_character_
    strand_downstream <- NA_character_

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
