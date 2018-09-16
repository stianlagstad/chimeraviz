#' Import results from a PRADA run into a list of Fusion objects.
#'
#' A function that imports the results from a PRADA run into a list of Fusion
#' objects.
#'
#' @param filename Filename for the PRADA results file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' pradaData <- system.file(
#'   "extdata",
#'   "PRADA.acc.fusion.fq.TAF.tsv",
#'   package = "chimeraviz")
#' fusions <- import_prada(pradaData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
import_prada <- function (filename, genome_version, limit) {

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
        "sample" = "character",
        "A_strand" = "integer",
        "B_strand" = "integer",
        "Discordant_n" = "integer",
        "JSR_n" = "character",
        "Junction" = "character"
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
      #message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      # End Exclude Linting
    }
  )

  # Set variables
  id                   <- NA
  inframe              <- NA
  fusion_tool          <- "prada"
  spanning_reads_count <- NA
  split_reads_count    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import prada-specific fields
    fusion_tool_specific_data <- list()

    # id for this fusion
    id <- as.character(i)

    # The format of report[[i, "Junction"]] should be:
    # "RBM14:11:66386032_BBS1:11:66283011,6". Split into two first:
    junction_data_gene_upstream <-
      strsplit(report[[i, "Junction"]], split = "_")[[1]][1]
    junction_data_gene_downstream <-
      strsplit(report[[i, "Junction"]], split = "_")[[1]][2]
    # The format of junctionDataGeneB should now be: "BBS1:11:66283011,6". So
    # remove everything after comma
    junction_data_gene_downstream <- substr(
      junction_data_gene_downstream,
      1,
      unlist(gregexpr(pattern = ",", junction_data_gene_downstream)) - 1)
    # Now split both on :
    junction_data_gene_upstream <-
      unlist(strsplit(junction_data_gene_upstream, split = ":"))
    junction_data_gene_downstream <-
      unlist(strsplit(junction_data_gene_downstream, split = ":"))

    # Chromosome names
    chromosome_upstream <-
      paste("chr", junction_data_gene_upstream[[2]], sep = "")
    chromosome_downstream <-
      paste("chr", junction_data_gene_downstream[[2]], sep = "")

    # Breakpoints
    breakpoint_upstream <- as.numeric(junction_data_gene_upstream[[3]])
    breakpoint_downstream <- as.numeric(junction_data_gene_downstream[[3]])

    # Strand
    strand_upstream <- if (report[[i, "A_strand"]] == 1) "+" else "-"
    strand_downstream <- if (report[[i, "B_strand"]] == 1) "+" else "-"

    # Number of supporting reads
    split_reads_count <- as.numeric(report[[i, "JSR_n"]])
    spanning_reads_count <- as.numeric(report[[i, "Discordant_n"]])

    # PRADA doesn't provide the fusion sequence
    junction_sequence_upstream <- Biostrings::DNAString()
    junction_sequence_downstream <- Biostrings::DNAString()

    # Gene names
    name_upstream <- junction_data_gene_upstream[[1]]
    name_downstream <- junction_data_gene_downstream[[1]]

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
