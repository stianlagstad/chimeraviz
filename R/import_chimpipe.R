#' Import results from a ChimPipe run into a list of Fusion objects.
#'
#' A function that imports the results from a ChimPipe run, typically from a
#' chimericJunctions_.txt file, into a list of Fusion objects.
#'
#' @param filename Filename for the ChimPipe results.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' chimpipefile <- system.file(
#'   "extdata",
#'   "chimericJunctions_MCF-7.txt",
#'   package="chimeraviz")
#' fusions <- import_chimpipe(chimpipefile, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.

#' @importFrom data.table fread
#'
#' @export
import_chimpipe <- function(filename, genome_version, limit) {

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
        "juncCoord" = "character",
        "gnTypesA" = "character",
        "gnTypesB" = "character",
        "gnNamesA" = "character",
        "gnNamesB" = "character",
        "nbSpanningReads" = "integer",
        "nbConsistentPE" = "integer"
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
  fusion_tool          <- "ChimPipe"
  spanning_reads_count <- NA
  split_reads_count    <- NA
  junction_sequence    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import defuse-specific fields
    fusion_tool_specific_data <- list() #TODO

    # Cluster id
    id <- as.character(i)

    # Is the downstream fusion partner in-frame?
    inframe <- NA

    parsed_junccord <- .chimpipe_parse_junccoord(report[[i, "juncCoord"]])

    # Strand
    strand_upstream <- parsed_junccord[3]
    strand_downstream <- parsed_junccord[6]

    # Number of supporting reads
    split_reads_count <- report[[i, "nbConsistentPE"]]
    spanning_reads_count <- report[[i, "nbSpanningReads"]]

    # No junction sequence is given from chimpipe
    junction_sequence_upstream <-
      Biostrings::DNAString()
    junction_sequence_downstream <-
      Biostrings::DNAString()

    # Breakpoints
    breakpoint_upstream <- tryCatch(
      as.numeric(parsed_junccord[2]),
      warning = function(w) {
        -1
      }
    )
    breakpoint_downstream <- tryCatch(
      as.numeric(parsed_junccord[5]),
      warning = function(w) {
        -1
      }
    )

    # Chromosome names
    chromosome_upstream <-
      parsed_junccord[1]
    chromosome_downstream <-
      parsed_junccord[4]

    # Gene names
    name_upstream <- report[[i, "gnNamesA"]]
    name_downstream <- report[[i, "gnNamesB"]]

    # Ensembl ids
    ensembl_id_upstream <- report[[i, "gnIdsA"]]
    ensembl_id_downstream <- report[[i, "gnIdsB"]]

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

.chimpipe_parse_junccoord <- function(junccoord) {
  if (class(junccoord) != "character") {
    stop("'junccoord' must be a character")
  }

  # Example input: chr12_9392928_+:chr11_75115716_+

  parts <- unlist(strsplit(junccoord, split = ":"))

  if (length(parts) != 2) {
    stop(paste0("Unable to parse 'junccoord':", junccoord))
  }

  parts_gene_upstream <- parts[1]
  parts_gene_downstream <- parts[2]

  upstream <- unlist(strsplit(parts_gene_upstream, split = "_"))
  downstream <- unlist(strsplit(parts_gene_downstream, split = "_"))
  
  if (length(upstream) != 3) {
    stop(paste0("Unable to parse 'junccoord':", junccoord))
  }
  if (length(downstream) != 3) {
    stop(paste0("Unable to parse 'junccoord':", junccoord))
  }

  c(upstream, downstream)
}
