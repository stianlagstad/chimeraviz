#' Import results from a Fusioncatcher run into a list of Fusion objects.
#'
#' A function that imports the results from a Fusioncatcher run, typically from
#' a final-list-candidate-fusion-genes.txt file, into a list of Fusion objects.
#'
#' @param filename Filename for the Fusioncatcher
#' final-list-candidate-fusion-genes.txt results file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' fusioncatcher833ke <- system.file(
#'   "extdata",
#'   "fusioncatcher_833ke_final-list-candidate-fusion-genes.txt",
#'   package = "chimeraviz")
#' fusions <- import_fusioncatcher(fusioncatcher833ke, "hg38", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
import_fusioncatcher <- function (filename, genome_version, limit) {

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
        "Gene_1_symbol(5end_fusion_partner)" = "character",
        "Gene_2_symbol(3end_fusion_partner)" = "character",
        "Spanning_pairs" = "integer",
        "Spanning_unique_reads" = "integer",
        "Fusion_point_for_gene_1(5end_fusion_partner)" = "character",
        "Fusion_point_for_gene_2(3end_fusion_partner)" = "character",
        "Gene_1_id(5end_fusion_partner)" = "character",
        "Gene_2_id(3end_fusion_partner)" = "character",
        "Fusion_sequence" = "character",
        "Predicted_effect" = "character"
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
  fusion_tool          <- "fusioncatcher"
  spanning_reads_count <- NA
  split_reads_count    <- NA
  junction_sequence    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import fusioncatcher-specific fields
    fusion_tool_specific_data <- list()

    # Set id for this fusion. Fusioncatcher does not have any unique identifier
    # for each row of the results file, so just use our counter.
    id <- as.character(i)

    # Is the downstream fusion partner in-frame?
    if (report[[i, "Predicted_effect"]] == "in-frame") {
      inframe <- TRUE
    } else {
      inframe <- FALSE
    }

    gene_upstream_position <-
      unlist(
        strsplit(
          report[[i, "Fusion_point_for_gene_1(5end_fusion_partner)"]],
          split = ":"
        )
      )
    gene_downstream_position <-
      unlist(
        strsplit(
          report[[i, "Fusion_point_for_gene_2(3end_fusion_partner)"]],
          split = ":"
        )
      )

    # Chromosome names
    chromosome_upstream <- paste("chr", gene_upstream_position[1], sep = "")
    chromosome_downstream <- paste("chr", gene_downstream_position[1], sep = "")

    # Breakpoints
    breakpoint_upstream <- as.numeric(gene_upstream_position[2])
    breakpoint_downstream <- as.numeric(gene_downstream_position[2])

    # Strand
    strand_upstream <- gene_upstream_position[3]
    strand_downstream <- gene_downstream_position[3]

    # Number of supporting reads
    split_reads_count <- report[[i, "Spanning_pairs"]]
    spanning_reads_count <- report[[i, "Spanning_unique_reads"]]

    # Get the fusion sequence. Split it into the part from gene1 and gene2
    junction_sequence <- strsplit(report[[i, "Fusion_sequence"]], split = "\\*")
    junction_sequence_upstream <-
      Biostrings::DNAString(junction_sequence[[1]][1])
    junction_sequence_downstream <-
      Biostrings::DNAString(junction_sequence[[1]][2])

    # Gene names
    name_upstream <- report[[i, "Gene_1_symbol(5end_fusion_partner)"]]
    name_downstream <- report[[i, "Gene_2_symbol(3end_fusion_partner)"]]

    # Ensembl ids
    ensembl_id_upstream <- report[[i, "Gene_1_id(5end_fusion_partner)"]]
    ensembl_id_downstream <- report[[i, "Gene_2_id(3end_fusion_partner)"]]

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
