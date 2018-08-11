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
      col_types_oncofuse <- readr::cols_only(
	"SAMPLE_ID" = col_skip(),
	"FUSION_ID" = readr::col_character(),
	"TISSUE" = col_skip(),
	"SPANNING_READS" = readr::col_integer(),
	"ENCOMPASSING_READS" = readr::col_integer(),
	"GENOMIC" = readr::col_character(),
	"5_FPG_GENE_NAME" = readr::col_character(),
	"5_IN_CDS" = col_skip(),
	"5_SEGMENT_TYPE" = col_skip(),
	"5_SEGMENT_ID" = col_skip(),
	"5_COORD_IN_SEGMENT" = col_skip(),
	"5_FULL_AA" = col_skip(),
	"5_FRAME" = col_skip(),
	"3_FPG_GENE_NAME" = readr::col_character(),
	"3_IN_CDS" = col_skip(),
	"3_SEGMENT_TYPE" = col_skip(),
	"3_SEGMENT_ID" = col_skip(),
	"3_COORD_IN_SEGMENT" = col_skip(),
	"3_FULL_AA" = col_skip(),
	"3_FRAME" = col_skip(),
	"FPG_FRAME_DIFFERENCE" = readr::col_integer(),
	"P_VAL_CORR" = col_skip(),
	"DRIVER_PROB" = readr::col_double(),
	"EXPRESSION_GAIN" = col_skip(),
	"5_DOMAINS_RETAINED" = col_skip(),
	"5_DOMAINS_RETAINED" = col_skip(),
	"5_DOMAINS_BROKEN" = col_skip(),
	"3_DOMAINS_BROKEN" = col_skip(),
	"5_PII_RETAINED" = col_skip(),
	"3_PII_RETAINED" = col_skip(),
	"CTF" = col_skip(),
	"G" = col_skip(),
	"H" = col_skip(),
	"K" = col_skip(),
	"P" = col_skip(),
	"TF" = col_skip()
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_oncofuse
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_oncofuse,
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
  fusion_tool          <- "oncofuse"
  spanning_reads_count <- NA
  split_reads_count    <- NA
  junction_sequence    <- NA
  strand_upstream      <- "+"
  strand_downstream    <- "+"

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

#    # Import oncofuse-specific fields
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

    # Strand
    #strand_upstream <- report[[i, "gene_strand1"]]
    #strand_downstream <- report[[i, "gene_strand2"]]

    # Number of supporting reads
    split_reads_count <- report[[i, "ENCOMPASSING_READS"]]
    spanning_reads_count <- report[[i, "SPANNING_READS"]]


    # Get the fusion sequence. Split it into the part from gene1 and gene2
    #junction_sequence <- strsplit(report[[i, "splitr_sequence"]], "\\|")
    #junction_sequence_upstream <-
    #  Biostrings::DNAString(junction_sequence[[1]][1])
    #junction_sequence_downstream <-
    #  Biostrings::DNAString(junction_sequence[[1]][2])

    # Breakpoints
    genomics <- report[[i, "GENOMIC"]]
    genomicVals <- unlist(strsplit(genomics, ">"))
    chromosome_upstream <- unlist(strsplit(genomicVals[1], ":"))[1]
    breakpoint_upstream <- as.numeric(unlist(strsplit(genomicVals[1], ":"))[2])
    chromosome_downstream <- unlist(strsplit(genomicVals[2], ":"))[1]
    breakpoint_downstream <- as.numeric(unlist(strsplit(genomicVals[2], ":"))[2])

    # Gene names
    name_upstream <- report[[i, "5_FPG_GENE_NAME"]]
    name_downstream <- report[[i, "3_FPG_GENE_NAME"]]

    # Ensembl REST API
    # Ensembl ids

    library(httr)
    library(jsonlite)
    library(xml2)
    server <- "https://grch37.rest.ensembl.org"
    #ext <- paste("/lookup/symbol/homo_sapiens/", name_upstream, "?format=condensed;db_type=core", sep="")
    ext <- paste("/lookup/symbol/homo_sapiens/", name_upstream, "?", sep="")
    r <- GET(paste(server, ext, sep=""), content_type("application/json"))
    stop_for_status(r)
    #print(content(r))

    ensemblObj <- content(r)
    ensembl_id_upstream <- ensemblObj$id
    
    ext <- paste("/sequence/region/human/", chromosome_upstream, ":", ensemblObj$start, "..", breakpoint_upstream, "?", sep="")
    r <- GET(paste(server, ext, sep=""), content_type("text/plain"))
    stop_for_status(r)
    #print(content(r))
    #junction_sequence_upstream <- content(r)
    junction_sequence_upstream <- Biostrings::DNAString(content(r))

    #ext <- paste("/lookup/symbol/homo_sapiens/", name_downstream, "?format=condensed;db_type=core", sep="")
    ext <- paste("/lookup/symbol/homo_sapiens/", name_downstream, "?", sep="")
    r <- GET(paste(server, ext, sep=""), content_type("application/json"))
    stop_for_status(r)
    ensemblObj <- content(r)
    ensembl_id_downstream <- ensemblObj$id

    ext <- paste("/sequence/region/human/", chromosome_downstream, ":", ensemblObj$start, "..", breakpoint_downstream, "?", sep="")
    r <- GET(paste(server, ext, sep=""), content_type("text/plain"))
    stop_for_status(r)
    #print(content(r))
    #junction_sequence_downstream <- content(r)
    junction_sequence_downstream <- Biostrings::DNAString(content(r))

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
