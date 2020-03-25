#' Import results from an Aeron run into a list of Fusion objects.
#'
#' A function that imports the results from an Aeron run into a list of Fusion
#' objects.
#' 
#' Note that the strands and breakpoint positions are not included in the
#' result files from Aeron. These have to be retrieved manually, using the
#' ensembl identifiers (which are included in the result files, and will be
#' available in the Fusion objects after importing).
#'
#' @param filename_fusion_support Filename for the Aeron result file
#' fusion_support..txt.
#' @param filename_fusion_transcript Filename for the Aeron result file
#' fusion_transcripts..txt.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' aeronfusionsupportfile <- system.file(
#'   "extdata",
#'   "aeron_fusion_support.txt",
#'   package="chimeraviz")
#' aeronfusiontranscriptfile <- system.file(
#'   "extdata",
#'   "aeron_fusion_transcripts.fa",
#'   package="chimeraviz")
#' fusions <- import_aeron(
#'   aeronfusionsupportfile,
#'   aeronfusiontranscriptfile,
#'   "hg19",
#'   3)
#' # This should import a list of 3 fusions described in Fusion objects.

#' @importFrom data.table fread
#'
#' @export
import_aeron <- function(
    filename_fusion_support,
    filename_fusion_transcript,
    genome_version,
    limit
  ) {

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

  # Try to read the fusion support report
  report <- withCallingHandlers({
      col_types <- c(
        "character",
        "integer"
      )
      if (missing(limit)) {
        # Read all lines
        data.table::fread(
          input = filename_fusion_support,
          colClasses = col_types,
          showProgress = FALSE,
          header = FALSE
        )
      } else {
          # Only read up to the limit
          data.table::fread(
            input = filename_fusion_support,
            colClasses = col_types,
            showProgress = FALSE,
            nrows = limit,
            header = FALSE
          )
      }
    },
    error = function(cond) {
      message(paste0("Reading ", filename_fusion_support, " caused an error: ", cond[[1]]))
      stop(cond)
    },
    warning = function(cond) {
      # Begin Exclude Linting
      message(paste0("Reading ", filename_fusion_support, " caused a warning: ", cond[[1]]))
      # End Exclude Linting
    }
  )

  # Read and parse the transcripts file
  fusion_transcript_map <-
    .import_aeron_transcript_data(filename_fusion_transcript)

  # Set variables
  id                   <- NA
  inframe              <- NA
  fusion_tool          <- "Aeron"
  spanning_reads_count <- NA
  split_reads_count    <- NA
  junction_sequence    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import defuse-specific fields
    fusion_tool_specific_data <- list()
    
    parsed_data <- .aeron_parse_data(report[[i, 1]])

    # Fusion id
    id <- parsed_data[1]

    # Get the data from the transcripts file
    transcript_data <- get(as.character(id), fusion_transcript_map)

    # Is the downstream fusion partner in-frame?
    inframe <- NA

    # Strand
    strand_upstream <- NA_character_
    strand_downstream <- NA_character_

    # Number of supporting reads
    split_reads_count <- NA_integer_
    spanning_reads_count <- tryCatch(
      as.numeric(parsed_data[4]),
      warning = function(w) {
        NA_integer_
      }
    )

    #Junction sequence
    junction_sequence_upstream <-
      Biostrings::DNAString(transcript_data[2])
    junction_sequence_downstream <-
      Biostrings::DNAString(transcript_data[3])

    # Breakpoints
    breakpoint_upstream <- -1
    breakpoint_downstream <- -1

    # Chromosome names
    chromosome_upstream <- NA_character_
    chromosome_downstream <- NA_character_

    # Gene names
    name_upstream <- NA_character_
    name_downstream <- NA_character_

    # Ensembl ids
    ensembl_id_upstream <- parsed_data[2]
    ensembl_id_downstream <- parsed_data[3]

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

.aeron_parse_data <- function(unparsed) {
  if (class(unparsed) != "character") {
    stop("'unparsed' must be a character")
  }

  # Example input: fusion_115_ENSG00000212907.2_291bp_ENSG00000198886.2_1340bp_50reads

  parts <- unlist(strsplit(unparsed, split = "_"))

  if (length(parts) != 7) {
    stop(paste0("Unable to parse 'unparsed':", unparsed))
  }

  fusion_id <- parts[2]
  upstream_id <- parts[3]
  downstream_id <- parts[5]
  spanning_reads <- gsub("[^0-9.-]", "", parts[7])

  c(fusion_id, upstream_id, downstream_id, spanning_reads)
}

.import_aeron_transcript_data <- function(
  filename_fusion_transcript
) {
  # Hashmap to hold all transcript data
  fusion_transcript_data_map <- new.env(hash=T, parent=emptyenv())

  # Read the .fa file
  dna = Biostrings::readDNAStringSet(filename_fusion_transcript)

  # Traverse the DNAStringSet
  for (i in 1:length(dna)) {
    entry <- dna[i]

    # Get the name
    fusion_name <- names(entry)
    # Split the fusion_name to get the fusion id
    fusion_name_parts <- unlist(strsplit(fusion_name, split = "_"))
    fusion_id <- fusion_name_parts[2]

    # Get the actual sequence
    actual_sequence <- toString(entry)

    # Split each on the "N" nucleotide
    actual_sequence_parts <- unlist(strsplit(actual_sequence, split = "N"))
    upstream_sequence <- actual_sequence_parts[1]
    downstream_sequence <- actual_sequence_parts[2]

    # Add this entry to our hashmap:
    assign(
      as.character(fusion_id),
      c(fusion_id, upstream_sequence, downstream_sequence),
      fusion_transcript_data_map
    )
  }

  fusion_transcript_data_map
}
