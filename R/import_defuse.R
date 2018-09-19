#' Import results from a deFuse run into a list of Fusion objects.
#'
#' A function that imports the results from a deFuse run, typically from a
#' results.filtered.tsv file, into a list of Fusion objects.
#'
#' @param filename Filename for the deFuse results .tsv file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#' 
#' @importFrom data.table fread
#'
#' @export
import_defuse <- function (filename, genome_version, limit) {

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
      col_types = c(
        "cluster_id" = "character",
        "splitr_sequence" = "character",
        "splitr_count" = "integer",
        "gene1" = "character",
        "gene2" = "character",
        "gene_chromosome1" = "character",
        "gene_chromosome2" = "character",
        "gene_name1" = "character",
        "gene_name2" = "character",
        "gene_strand1" = "character",
        "gene_strand2" = "character",
        "genomic_break_pos1" = "integer",
        "genomic_break_pos2" = "integer",
        "span_count" = "integer",
        "orf" = "character",
        "probability" = "numeric"
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
  fusion_tool          <- "defuse"
  spanning_reads_count <- NA
  split_reads_count    <- NA
  junction_sequence    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import defuse-specific fields
    fusion_tool_specific_data <- list()
    fusion_tool_specific_data[["probability"]] <- report[[i, "probability"]]

    # Cluster id
    id <- report[[i, "cluster_id"]]

    # Is the downstream fusion partner in-frame?
    if (report[[i, "orf"]] == "Y") {
      inframe <- TRUE
    } else {
      inframe <- FALSE
    }

    # Strand
    strand_upstream <- report[[i, "gene_strand1"]]
    strand_downstream <- report[[i, "gene_strand2"]]

    # Number of supporting reads
    split_reads_count <- report[[i, "splitr_count"]]
    spanning_reads_count <- report[[i, "span_count"]]

    # Get the fusion sequence. Split it into the part from gene1 and gene2
    junction_sequence <- strsplit(report[[i, "splitr_sequence"]], "\\|")
    junction_sequence_upstream <-
      Biostrings::DNAString(junction_sequence[[1]][1])
    junction_sequence_downstream <-
      Biostrings::DNAString(junction_sequence[[1]][2])

    # Breakpoints
    breakpoint_upstream <- report[[i, "genomic_break_pos1"]]
    breakpoint_downstream <- report[[i, "genomic_break_pos2"]]

    # Chromosome names
    chromosome_upstream <-
      paste("chr", report[[i, "gene_chromosome1"]], sep = "")
    chromosome_downstream <-
      paste("chr", report[[i, "gene_chromosome2"]], sep = "")

    # Gene names
    name_upstream <- report[[i, "gene_name1"]]
    name_downstream <- report[[i, "gene_name2"]]

    # Ensembl ids
    ensembl_id_upstream <- report[[i, "gene1"]]
    ensembl_id_downstream <- report[[i, "gene2"]]

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
