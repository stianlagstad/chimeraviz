#' Fetch reads from fastq files
#'
#' This function will fetch read sequences from fastq files and put them into
#' new fastq files.
#'
#' Note: This function runs (read only) bash commands on your system. Therefore
#' the function will only work on a unix system.
#'
#' @param reads List of read IDs that is to be fetched.
#'
#' @param fastq_file_in1 First fastq file to search in.
#' @param fastq_file_in2 Second fastq file to seach in.
#'
#' @param fastq_file_out1 First fastq file with results.
#' @param fastq_file_out2 Second fastq file with results.
#'
#' @return The files fastqFileOut1 and fastqFileOut2 populated with the
#' specified reads.
#'
#' @examples
#' \dontrun{
#' # fastq files that has the supporting reads
#' fastq1 <- system.file("extdata", "reads.1.fq", package="chimeraviz")
#' fastq2 <- system.file("extdata", "reads.2.fq", package="chimeraviz")
#' # Which read ids to extract
#' reads <- c(
#'   "13422259", "19375605", "29755061",
#'   "31632876", "32141428", "33857245")
#' # Extract the actual reads and put them in the tmp files "fastqFileOut1" and
#' # "fastqFileOut2"
#' fastqFileOut1 <- tempfile(pattern = "fq1", tmpdir = tempdir())
#' fastqFileOut2 <- tempfile(pattern = "fq2", tmpdir = tempdir())
#' fetch_reads_from_fastq(reads, fastq1, fastq2,
#'     fastqFileOut1, fastqFileOut2)
#' # We now have the reads supporting fusion 5267 in the two files.
#' }
#'
#' @export
fetch_reads_from_fastq <- function(
  reads,
  fastq_file_in1,
  fastq_file_in2,
  fastq_file_out1,
  fastq_file_out2
  ) {

  # Since this function use the bash function egrep to extract reads, give a
  # warning if we're running on windows
  if (Sys.info()["sysname"] == "Windows") {
    stop(paste("This function uses the bash function egrep. It looks like",
               "you're running on windows, so this function will terminate."))
  }

  if (is.vector(reads, mode = "character") == FALSE) {
    stop("reads should be a character vector of read ids")
  }

  if (file.exists(fastq_file_in1) == FALSE ||
      file.exists(fastq_file_in2) == FALSE) {
    stop("Invalid fastq input files")
  }

  # The command below will extract the sequences for ids 11, 22, and 33 from the
  # fastq file reads.fastq

  # $ egrep -A 3 '@11|22|33' reads.fastq | sed '/^--$/d'

  # By default, the -A parameter inserts a "--" between matches, so the sed part
  # above removes that

  # First build the query with the read ids
  query <- paste("@",
                 paste(reads, collapse = "|"),
                 sep = "")

  # Then create commands for each fastq file
  command1 <- paste("egrep -A 3 '",
                    query,
                    "' ",
                    shQuote(fastq_file_in1),
                    " | sed '/^--$/d'", # Exclude Linting
                    sep = "")
  command2 <- paste("egrep -A 3 '",
                    query,
                    "' ",
                    shQuote(fastq_file_in2),
                    " | sed '/^--$/d'", # Exclude Linting
                    sep = "")

  # We now have two commands:
  # egrep -A 1 '@11|22|33' reads.1.fastq | sed '/^--$/d'
  # egrep -A 1 '@11|22|33' reads.2.fastq | sed '/^--$/d'

  # Run the command and capture output
  supporting_reads_fq1 <- system(command1, intern = TRUE)
  supporting_reads_fq2 <- system(command2, intern = TRUE)

  # Write these to each their own files
  write(supporting_reads_fq1, file = fastq_file_out1)
  write(supporting_reads_fq2, file = fastq_file_out2)

  # Check file sizes
  if (file.info(fastq_file_out1)$size <= 10) {
    warning(paste("It looks like we couldn't find any of the reads in the",
                  "fastq file. Check your arguments."))
  }
}

#' Write fusion junction sequence to a fasta file
#'
#' This function will write the fusion sequence to a fasta file, using
#' Biostring::writeXStringSet() .
#'
#' @param fusion The Fusion object we want to create a fasta file from.
#' @param filename The filename to write to.
#'
#' @return Writes the fusion junction sequence to the given filename.
#'
#' @examples
#' # Import the filtered defuse results
#' defuse833keFiltered <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833keFiltered, "hg19", 1)
#' # Get a specific fusion
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Create temporary file to hold the fusion sequence
#' fastaFileOut <- tempfile(pattern = "fusionSequence", tmpdir = tempdir())
#' # Write fusion sequence to file
#' write_fusion_reference(fusion, fastaFileOut)
#'
#' @export
write_fusion_reference <- function(fusion, filename) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  # First put the fusion junction sequence in a DNAStringSet object
  fusion_sequence <- Biostrings::DNAStringSet(
    x = c(fusion@gene_upstream@junction_sequence,
          fusion@gene_downstream@junction_sequence))

  # Give an error if the length of the fusionSequence is 0:
  if (nchar(fusion_sequence) == 0) {
    stop(
      paste0(
        "The fusion sequence length is zero, so the fusion reference sequence",
        " cannot be written."
      )
    )
  }

  # Set sequence name to chrNA, since this is a sequence created from a fusion
  # event (i.e. not a sequence from a real chromosome). The "chrNA" name will
  # make Gviz happy.
  names(fusion_sequence) <- "chrNA"
  # Write to file
  Biostrings::writeXStringSet(fusion_sequence, filename)
}

#' Get ensembl ids for a fusion object
#'
#' This function will get the ensembl ids from the org.Hs.eg.db/org.Mm.eg.db
#' package given the gene names of the fusion event.
#'
#' @param fusion The Fusion object we want to get ensembl ids for.
#'
#' @return The Fusion object with Ensembl ids set.
#'
#' @import org.Hs.eg.db org.Mm.eg.db
#'
#' @examples
#' # Import the filtered defuse results
#' defuse833keFiltered <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833keFiltered, "hg19", 1)
#' # Get a specific fusion
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # See the ensembl ids:
#' partner_gene_ensembl_id(upstream_partner_gene(fusion))
#' # [1] "ENSG00000180198"
#' partner_gene_ensembl_id(downstream_partner_gene(fusion))
#' # [1] "ENSG00000162639"
#' # Reset the fusion objects ensembl ids
#' partner_gene_ensembl_id(upstream_partner_gene(fusion)) <- ""
#' partner_gene_ensembl_id(downstream_partner_gene(fusion))  <- ""
#' # Get the ensembl ids
#' fusion <- get_ensembl_ids(fusion)
#' # See that we now have the same ensembl ids again:
#' partner_gene_ensembl_id(upstream_partner_gene(fusion))
#' # [1] "ENSG00000180198"
#' partner_gene_ensembl_id(downstream_partner_gene(fusion))
#' # [1] "ENSG00000162639"
#'
#' @export
get_ensembl_ids <- function(fusion) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  if (startsWith(fusion@genome_version, "hg")) {
    annotation_db <- org.Hs.eg.db # Exclude Linting
  } else if (startsWith(fusion@genome_version, "mm")) {
    annotation_db <- org.Mm.eg.db # Exclude Linting
  } else {
    stop("Unsupported genome version")
  }

  result <- AnnotationDbi::select(
    annotation_db,
    keys = c(fusion@gene_upstream@name, fusion@gene_downstream@name),
    keytype = "ALIAS",
    columns = c("ALIAS", "ENSEMBL"))

  upstream_gene_result <-
    result[which(result$ALIAS == fusion@gene_upstream@name), ]
  downstream_gene_result <-
    result[which(result$ALIAS == fusion@gene_downstream@name), ]

  # Stop execution if no results
  if (any(is.na(upstream_gene_result$ENSEMBL))) {
    stop(
      paste(
        "Could not find Ensembl id for ",
        fusion@gene_upstream@name,
        ". ",
        "If you know the id, add it manually with ",
        "fusion@gene_upstream@ensembl_id <- \"ensemblId\"",
        sep = ""
      )
    )
  }
  if (any(is.na(downstream_gene_result$ENSEMBL))) {
    stop(
      paste(
        "Could not find Ensembl id for ",
        fusion@gene_downstream@name,
        ". ",
        "If you know the id, add it manually with ",
        "fusion@gene_downstream@ensembl_id <- \"ensemblId\"",
        sep = ""
      )
    )
  }

  # Store ensembl ids in fusion object
  fusion@gene_upstream@ensembl_id <- upstream_gene_result[, 2]
  fusion@gene_downstream@ensembl_id <- downstream_gene_result[, 2]

  # Return updated fusion object
  fusion
}

#' Split GRanges object based on cds
#'
#' This function will look for ranges (exons) in the GRanges object that has the
#' coding DNA sequence starting or stopping within it. If found, these exons are
#' split, and each exon in the GRanges object will be tagged as either
#' "protein_coding", "5utr", or "3utr". The returned GRanges object will have
#' feature values set in mcols(gr)$feature reflecting this.
#'
#' @param gr The GRanges object we want to split and tag with feature info.
#'
#' @return An updated GRanges object with feature values set.
#'
#' @examples
#' # Load fusion data and choose a fusion object:
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Create edb object
#' edbSqliteFile <- system.file(
#'   "extdata",
#'   "Homo_sapiens.GRCh37.74.sqlite",
#'   package="chimeraviz")
#' edb <- ensembldb::EnsDb(edbSqliteFile)
#' # Get all exons for all transcripts in the genes in the fusion transcript
#' allTranscripts <- ensembldb::exonsBy(
#'   edb,
#'   filter = list(
#'     AnnotationFilter::GeneIdFilter(
#'       list(
#'         partner_gene_ensembl_id(upstream_partner_gene(fusion)),
#'         partner_gene_ensembl_id(downstream_partner_gene(fusion))))),
#'   columns = c(
#'     "gene_id",
#'     "gene_name",
#'     "tx_id",
#'     "tx_cds_seq_start",
#'     "tx_cds_seq_end",
#'     "exon_id"))
#' # Extract one of the GRanges objects
#' gr <- allTranscripts[[1]]
#' # Check how many ranges there are here
#' length(gr)
#' # Should be 9 ranges
#' # Split the ranges containing the cds start/stop positions and add feature
#' # values:
#' gr <- split_on_utr_and_add_feature(gr)
#' # Check the length again
#' length(gr)
#' # Should be 11 now, as the range containing the cds_strat position and the
#' # range containing the cds_stop position has been split into separate ranges
#'
#' @importFrom S4Vectors mcols
#'
#' @export
split_on_utr_and_add_feature <- function(gr) {

  # Check if we got a valid GRanges object
  if (class(gr) != "GRanges") {
    stop("gr argument must be an object of type GRanges")
  }

  # Check that we have the tx_cds_seq_start and tx_cds_seq_end values


  # check that gr has a tx_cds_seq_start and tx_cds_seq_end mcols value

  # Get cds start and cds end
  cds_start <- S4Vectors::mcols(gr)$tx_cds_seq_start[[1]]
  cds_end <- S4Vectors::mcols(gr)$tx_cds_seq_end[[1]]

  # find exon containing start
  cds_start_exon <- cds_start > start(gr) & cds_start < end(gr)

  if (any(cds_start_exon)) {
    # Create two copies of the exon in which the cds starts
    first <- gr[cds_start_exon]
    second <- gr[cds_start_exon]
    # Update ranges for the new objects
    end(first) <- cds_start - 1
    start(second) <- cds_start
    # Remove the original range
    gr <- gr[!cds_start_exon]
    # Add new ranges
    gr <- append(gr, first)
    gr <- append(gr, second)
    # Sort again
    gr <- sort(gr)
  }

  # find exon containing end
  cds_end_exon <- cds_end > start(gr) & cds_end < end(gr)

  if (any(cds_end_exon)) {
    # Create two copies of the exon in which the cds ends
    first <- gr[cds_end_exon]
    second <- gr[cds_end_exon]
    # Update ranges for the new objects
    end(first) <- cds_end
    start(second) <- cds_end + 1
    # Remove the original range
    gr <- gr[!cds_end_exon]
    # Add new ranges
    gr <- append(gr, first)
    gr <- append(gr, second)
    # Sort again
    gr <- sort(gr)
  }

  # add features
  S4Vectors::mcols(gr)$feature <- "protein_coding"
  # 5utr:
  #   plus strand:
  #     end(gr) < cds_start # Exclude Linting
  #   minus strand:
  #     start(gr) > cds_end # Exclude Linting
  if (length(gr[as.character(strand(gr)) == "+" & end(gr) < cds_start |
                as.character(strand(gr)) == "-" & start(gr) > cds_end])) {
    S4Vectors::mcols(
      gr[as.character(strand(gr)) == "+" & end(gr) < cds_start |
           as.character(strand(gr)) == "-" & start(gr) > cds_end]
      )$feature <- "5utr"
  }
  # 3utr:
  #   plus strand:
  #     start(gr) > cds_end # Exclude Linting
  #   minus strand:
  #     end(gr) < cds_start # Exclude Linting
  if (length(gr[as.character(strand(gr)) == "+" & start(gr) > cds_end |
                as.character(strand(gr)) == "-" & end(gr) < cds_start])) {
    S4Vectors::mcols(
      gr[as.character(strand(gr)) == "+" & start(gr) > cds_end |
           as.character(strand(gr)) == "-" & end(gr) < cds_start]
      )$feature <- "3utr"
  }

  gr
}

#' Retrieves transcripts for partner genes in a Fusion object using Ensembldb
#'
#' This function will check where in the transcript (the GRanges object) the
#' fusion breakpoint is located, and return either "exonBoundary", "withinExon",
#' "withinIntron", or "intergenic".
#'
#' @param gr The GRanges object containing the transcript to be checked.
#' @param fusion The fusion object used to check the transcript.
#'
#' @return Either "exonBoundary", "withinExon", "withinIntron", or "intergenic"
#' depending on where in the transcript the breakpoint hits.
#'
#' @examples
#' # Load fusion data and choose a fusion object:
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Create edb object
#' edbSqliteFile <- system.file(
#'   "extdata",
#'   "Homo_sapiens.GRCh37.74.sqlite",
#'   package="chimeraviz")
#' edb <- ensembldb::EnsDb(edbSqliteFile)
#' # Get all exons for all transcripts in the genes in the fusion transcript
#' allTranscripts <- ensembldb::exonsBy(
#'   edb,
#'   filter = list(
#'     AnnotationFilter::GeneIdFilter(
#'       list(
#'         partner_gene_ensembl_id(upstream_partner_gene(fusion)),
#'         partner_gene_ensembl_id(downstream_partner_gene(fusion))))),
#'   columns = c(
#'     "gene_id",
#'     "gene_name",
#'     "tx_id",
#'     "tx_cds_seq_start",
#'     "tx_cds_seq_end",
#'     "exon_id"))
#' # Extract one of the GRanges objects
#' gr <- allTranscripts[[1]]
#' # Check where in the transcript the fusion breakpoint hits
#' decide_transcript_category(gr, fusion)
#' # "exonBoundary"
#' # Check another case
#' gr <- allTranscripts[[3]]
#' decide_transcript_category(gr, fusion)
#' # "withinIntron"
#'
#' @importFrom S4Vectors mcols
#'
#' @export
decide_transcript_category <- function(gr, fusion) {

  # Check if we got a valid GRanges object
  if (class(gr) != "GRanges") {
    stop("gr argument must be an object of type GRanges")
  }

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  # 4 possible cases:
  #
  # 1: geneA is on the + strand
  #    Exon boundary event if breakpointA equals an exon_end position
  # 2: geneA is on the - strand
  #    Exon boundary event if breakpointA equals an exon_start position
  # 3: geneB is on the + strand
  #    Exon boundary event if breakpointB equals an exon_start position
  # 4: geneB is on the - strand
  #    Exon boundary event if breakpointB equals an exon_end position

  # Helper variables
  upstream_gene <-
    S4Vectors::mcols(gr)$gene_id[[1]] == fusion@gene_upstream@ensembl_id
  breakpoint <-
    if (all(S4Vectors::mcols(gr)$gene_id == fusion@gene_upstream@ensembl_id)) {
      fusion@gene_upstream@breakpoint
    } else {
      fusion@gene_downstream@breakpoint
    }
  exon_start_positions <- GenomicRanges::start(gr)
  exon_end_positions <- GenomicRanges::end(gr)

  # Check whether or not the breakpoint is within the transcript
  if (all(breakpoint < exon_start_positions) ||
      all(exon_end_positions < breakpoint)) {
    return("intergenic")
  }

  # Assuming only single-stranded transcripts
  if (as.character(strand(gr)[1]) == "+") {

    # To decide whether or not the breakpoint occurs at an exon boundary, we
    # have to be vary careful about both whether this is the upstream/downstrem
    # fusion partner gene
    if (upstream_gene) {
      # GeneA
      if (any(exon_end_positions == breakpoint)) {
        # exon boundary
        return("exonBoundary")
      }
    } else {
      # GeneB
      if (any(exon_start_positions == breakpoint)) {
        # exon boundary
        return("exonBoundary")
      }
    }

  } else {

    # To decide whether or not the breakpoint occurs at an exon boundary, we
    # have to be vary careful about both whether this is the upstream/downstrem
    # fusion partner gene
    if (all(upstream_gene)) {
      # GeneA
      if (any(exon_start_positions == breakpoint)) {
        # exon boundary
        return("exonBoundary")
      }
    } else {
      # GeneB
      if (any(exon_end_positions == breakpoint)) {
        # exon boundary
        return("exonBoundary")
      }
    }

  }

  # Within exon if end < breakpoint & breakpoint < start
  if (any(start(gr) < breakpoint & breakpoint < end(gr))) {
    return("withinExon")
  } else {
    return("withinIntron")
  }
}

#' Retrieves transcripts for partner genes in a Fusion object using Ensembldb
#'
#' This function will retrieve transcripts for both genes in a fusion. It will
#' check all transcripts and decide for each transcript if the fusion breakpoint
#' happens at 1) an exon boundary, 2) within an exon, or 3) within an intron.
#' This is done because fusions happening at exon boundaries are more likely to
#' produce biologically interesting gene products. The function returns an
#' updated Fusion object, where the fusion@gene_upstream@transcriptsX slots are set with
#' transcript information.
#'
#' @param fusion The fusion object to find transcripts for.
#' @param edb The edb object used to fetch data from.
#'
#' @return An updated fusion object with transcript data stored.
#'
#' @examples
#' # Load fusion data and choose a fusion object:
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Create edb object
#' edbSqliteFile <- system.file(
#'   "extdata",
#'   "Homo_sapiens.GRCh37.74.sqlite",
#'   package="chimeraviz")
#' edb <- ensembldb::EnsDb(edbSqliteFile)
#' # Add transcripts data to fusion object
#' fusion <- get_transcripts_ensembl_db(fusion, edb)
#' # The transcripts are now accessible through fusion@gene_upstream@transcripts and
#' # fusion@gene_downstream@transcripts .
#'
#' @importFrom S4Vectors mcols
#' @importFrom ensembldb exonsBy
#' @importFrom AnnotationFilter GeneIdFilter
#'
#' @export
get_transcripts_ensembl_db <- function(fusion, edb) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  # Check if we got a valid edb
  if (class(edb) != "EnsDb") {
    stop("edb argument must be an object of type EnsDb")
  }

  # Fetch ensembl gene ids if we don't already have them
  if (
    is.na(fusion@gene_upstream@ensembl_id) ||
    is.na(fusion@gene_downstream@ensembl_id)
  ) {
    fusion <- get_ensembl_ids(fusion)
  }

  # Get all exon information
  all_transcripts <- ensembldb::exonsBy(
    edb,
    filter = list(
      AnnotationFilter::GeneIdFilter(
        list(
          fusion@gene_upstream@ensembl_id,
          fusion@gene_downstream@ensembl_id))),
    columns = c(
      "gene_id",
      "gene_name",
      "tx_id",
      "tx_cds_seq_start",
      "tx_cds_seq_end",
      "exon_id"))

  # Fail if no transcripts were found
  if (length(all_transcripts) == 0) {
    stop(paste(
      "No transcripts available for the genes ",
      fusion@gene_upstream@name,
      " and ",
      fusion@gene_downstream@name,
      ".",
      sep = ""))
  }

  # Go through each transcript in the GRangesList and
  grangeslist_upstream <- GRangesList()
  grangeslist_downstream <- GRangesList()
  for (i in 1:length(all_transcripts)) {

    # Extract the GRanges object
    gr <- all_transcripts[[i]]

    # Add $transcript as a metadata column
    S4Vectors::mcols(gr)$transcript <- names(all_transcripts[i])
    # Add $symbol as a metadata column
    S4Vectors::mcols(gr)$symbol <- names(all_transcripts[i])

    # Add utr/coding S4Vectors::mcols$feature data for each region
    if (is.na(S4Vectors::mcols(gr)$tx_cds_seq_start[[1]])) {
      # The transcript is not coding. Set all to "utr"
      S4Vectors::mcols(gr)$feature <- "utr"
    }else {
      # The transcript is coding. The function below will set
      # S4Vectors::mcols(gr)$feature = "utr5"/"protein_coding"/"utr3" for each
      # region
      gr <- split_on_utr_and_add_feature(gr)
    }

    # Create new GRangesList
    grl <- GRangesList(gr)

    # Tag each transcript with either "exonBoundary", "withinExon",
    # "withinIntron", or "intergenic", depending on where in the transcript the
    # fusion breakpoint hits
    S4Vectors::mcols(grl)$transcript_category <- decide_transcript_category(
      gr,
      fusion)

    # Set S4Vectors::mcols$transcript and name for the grl
    S4Vectors::mcols(grl)$transcript <- names(all_transcripts[i])
    names(grl) <- names(all_transcripts[i])

    # Add S4Vectors::mcols$coding to signify whether it is a coding transcript
    # at all
    S4Vectors::mcols(grl)$coding <-
      !all(is.na(S4Vectors::mcols(gr)$tx_cds_seq_start))

    # Append the GRanges object to the correct GRangesList
    if (S4Vectors::mcols(gr)$gene_id[[1]] == fusion@gene_upstream@ensembl_id) {
      grangeslist_upstream <- append(grangeslist_upstream, grl)
    } else {
      grangeslist_downstream <- append(grangeslist_downstream, grl)
    }
  }

  # Add transcripts to fusion object
  fusion@gene_upstream@transcripts <- grangeslist_upstream
  fusion@gene_downstream@transcripts <- grangeslist_downstream

  # Warn if no transcripts were found for one of the genes
  if (length(fusion@gene_upstream@transcripts) == 0) {
    warning(paste(
      "No transcripts available for the upstream gene ",
      fusion@gene_upstream@name,
      " available.",
      sep = ""))
  }
  if (length(fusion@gene_upstream@transcripts) == 0) {
    warning(paste(
      "No transcripts available for the downstream gene ",
      fusion@gene_downstream@name,
      " available.",
      sep = ""))
  }

  # In case the fusion object doesn't have the strands set (as is the case for
  # JAFFA), set the strands now:
  if (
    fusion@gene_upstream@strand == "*" || fusion@gene_downstream@strand == "*"
  ) {
    fusion@gene_upstream@strand <- as.character(
      strand(fusion@gene_upstream@transcripts[[1]][1]))
    fusion@gene_downstream@strand <- as.character(
      strand(fusion@gene_downstream@transcripts[[1]][1]))
  }

  # Return fusion object
  fusion

}

#' Add fusion reads alignment to fusion object
#'
#' This function lets you add a fusion read alignment file to a fusion object.
#' If you've mapped the reads supporting a fusion against the fusion junction
#' sequence, and have the resulting bamfile, use this function to add the
#' information (as a Gviz::GAlignmentPairs object) to the fusion object.
#'
#' @param fusion The fusion object to add a genomic alignment to.
#' @param bamfile The bam file containing the fusion reads plotted to the fusion
#' sequence.
#'
#' @return An updated fusion object with fusion@fusion_reads_alignment set.
#'
#' @examples
#' # Load data
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' # Find the specific fusion we have aligned reads for
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Get reference to the bamfile with the alignment data
#' bamfile5267 <- system.file(
#'   "extdata",
#'   "5267readsAligned.bam",
#'   package="chimeraviz")
#' # Add the bam file of aligned fusion reads to the fusion object
#' fusion <- add_fusion_reads_alignment(fusion, bamfile5267)
#'
#' @export
add_fusion_reads_alignment <- function(fusion, bamfile) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  fusion_reads_alignment <- Gviz::AlignmentsTrack(
    bamfile,
    isPaired = TRUE,
    # Setting chromosome to chrNA because this is a fusion sequence not found in
    # any reference genome.
    chromosome = "chrNA",
    name = "Fusion Reads",
    genome = fusion@genome_version)

  # Return new fusion object, now with the fusion read alignment
  fusion@fusion_reads_alignment <- fusion_reads_alignment
  fusion
}

#' Coerce Fusion object to data.frame
#'
#' This function is used in create_fusion_report() to convert Fusion objects to a
#' data.frame-format.
#'
#' @param fusion The Fusion object to coerce.
#'
#' @return A data.frame with the fusion object.
#'
#' @seealso create_fusion_report
#'
#' @examples
#' # Load data
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' # Find the fusion object to create a data frame from
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Create the data frame
#' dfFusion <- fusion_to_data_frame(fusion)
#'
#' @export
fusion_to_data_frame <- function(fusion) {

  # Check if we got a list of fusion objects
  if (class(fusion) != "Fusion") {
    stop("fusions argument must be a list of Fusion objects")
  }

  df <- data.frame(
    fusion@id,
    fusion@gene_upstream@name,
    fusion@gene_upstream@ensembl_id,
    fusion@gene_upstream@breakpoint,
    fusion@gene_downstream@name,
    fusion@gene_downstream@ensembl_id,
    fusion@gene_downstream@breakpoint,
    fusion@split_reads_count,
    fusion@spanning_reads_count)
  names(df) <- c(
    "id",
    "gene_upstream",
    "ensembl_upstream",
    "breakpoint_upstrea,",
    "gene_downstream",
    "ensembl_downstream",
    "breakpoint_downstream",
    "Split Reads",
    "Spanning Reads")
  df
}

#' Select which transcript to use (for plotting) for a GenePartner object
#'
#' This function takes a GenePartner object and creates a transcript data.frame
#' with transcript information, including only the transcripts given by the
#' parameter which_transcripts
#'
#' select_transcript() selects which transcript to create by this prioritization:
#'
#' 1. Exon boundary transcripts.
#' 2. Within exon transcripts.
#' 3. Within intron transcripts.
#' 4. Intergenic transcripts.
#'
#' @param gene_partner The GenePartner object to select a transcript for.
#' @param which_transcripts This character vector decides which transcripts are
#' to be plotted. Can be "exonBoundary", "withinExon", "withinIntron",
#' "intergenic", or a character vector with specific transcript ids. Default
#' value is "exonBoundary".
#'
#' @return A data.frame with transcript data.
#'
#' @examples
#' # Load data and example fusion event
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # Load edb
#' edbSqliteFile <- system.file(
#'   "extdata",
#'   "Homo_sapiens.GRCh37.74.sqlite",
#'   package="chimeraviz")
#' edb <- ensembldb::EnsDb(edbSqliteFile)
#' # Get transcripts
#' fusion <- get_transcripts_ensembl_db(fusion, edb)
#' # Select transcript
#' transcriptsA <- select_transcript(upstream_partner_gene(fusion))
#'
#' @export
select_transcript <- function(
  gene_partner,
  which_transcripts = "exonBoundary") {

  # Check if we got a PartnerGene object
  if (class(gene_partner) != "PartnerGene") {
    stop("genePartner argument must be an object of type PartnerGene")
  }

  # Does the PartnerGene have any transcripts?
  if (isEmpty(gene_partner@transcripts)) {
    stop("genePartner has no transcripts. See get_transcripts_ensembl_db()")
  }

  # Either select from one of these four categories
  # 1. Exon boundary transcripts.
  # 2. Within exon transcripts.
  # 3. Within intron transcripts.
  # 4. Intergenic transcripts.
  #
  # or select the specific transcripts given by which_transcripts.

  transcript_categories <- c(
    "exonBoundary",
    "withinExon",
    "withinIntron",
    "intergenic")

  if (which_transcripts[[1]] %in% transcript_categories) {

    message(paste0("Selecting transcripts for ", gene_partner@name, ".."))

    # If the user has chosen one of the four transcript categories, then we
    # want to check whether or not such transcripts exist. If they exist,
    # simply return them. If they don't exist, go on to try the other
    # categories. Try the wanted category first:
    transcripts_of_wanted_category <-
      gene_partner@transcripts[
        mcols(gene_partner@transcripts)$
          transcript_category == which_transcripts[[1]]
      ]
    if (length(transcripts_of_wanted_category) > 0) {
      message(paste0("..found transcripts of type ", which_transcripts[[1]]))
      return(transcripts_of_wanted_category)
    }
    # Check the remaining categories
    remaining_categories <- transcript_categories[
      transcript_categories != which_transcripts[[1]]
    ]
    for (transcript_category in remaining_categories) {
      transcripts_in_this_category <- gene_partner@transcripts[
        mcols(gene_partner@transcripts)$
          transcript_category == transcript_category
      ]
      if (length(transcripts_in_this_category) > 0) {
        message(paste0("..found transcripts of type ", transcript_category))
        return(transcripts_in_this_category)
      }
    }
  }

  # At this point the user wants specific transcripts. Get the transcripts that
  # we have for this genePartner.

  specific_transcripts <- gene_partner@transcripts[
    names(gene_partner@transcripts) %in% which_transcripts]
  if (length(specific_transcripts) > 0) {
    return(specific_transcripts)
  }

  stop("The specific transcripts could not be found")
}

# Check that there's at least one transcript that has the fusion breakpoint
# within the transcript.
.check_that_breakpoints_are_within_transcripts <- function(
  fusion,
  transcripts_upstream,
  transcripts_downstream) {
  if (!any(start(transcripts_upstream) < fusion@gene_upstream@breakpoint) &
      any(fusion@gene_upstream@breakpoint < end(transcripts_upstream))) {
    stop(paste0(
      "None of the transcripts given for gene A has the fusion breakpoint ",
      "within them. This plot cannot be created with the given transcripts."))
  }
  if (!any(start(transcripts_downstream) < fusion@gene_downstream@breakpoint) &
      any(fusion@gene_downstream@breakpoint < end(transcripts_downstream))) {
    stop(paste0(
      "None of the transcripts given for gene B has the fusion breakpoint ",
      "within them. This plot cannot be created with the given transcripts."))
  }
}

# Check that the transcripts have a breakpoint exon
.check_that_transcripts_have_breakpoint_exons <- function(
  fusion,
  transcript_upstream,
  transcript_downstream) {

  if (fusion@gene_upstream@strand == "+") {
    if (
      length(
        transcript_upstream[
          end(ranges(transcript_upstream)) == fusion@gene_upstream@breakpoint
        ]
      ) == 0
    ) {
      stop(paste(
        "The transcript for gene A doesn't have an exon boundary matching the",
        "fusion breakpoint. The plot cannot be created."))
    }
  } else {
    if (
      length(
        transcript_upstream[
          start(ranges(transcript_upstream)) == fusion@gene_upstream@breakpoint
        ]
      ) == 0
    ) {
      stop(paste(
        "The transcript for gene A doesn't have an exon boundary matching the",
        "fusion breakpoint. The plot cannot be created."))
    }
  }
  if (fusion@gene_downstream@strand == "+") {
    if (
      length(
        transcript_downstream[
          start(ranges(transcript_downstream)) ==
          fusion@gene_downstream@breakpoint
        ]
      ) == 0
    ) {
      stop(paste(
        "The transcript for gene B doesn't have an exon boundary matching the",
        "fusion breakpoint. The plot cannot be created."))
    }
  } else {
    if (
      length(
        transcript_downstream[
          end(ranges(transcript_downstream)) ==
          fusion@gene_downstream@breakpoint
        ]
      ) == 0
    ) {
      stop(paste(
        "The transcript for gene B doesn't have an exon boundary matching the",
        "fusion breakpoint. The plot cannot be created."))
    }
  }
}

#' Remove introns and shift exons leftward
#'
#' This function takes a GRanges object and moves each IRanges object within
#' next to each other starting at 1. This effectively removes the introns from
#' the GRanges object.
#'
#' @param transcript The GRanges object to remove introns from.
#'
#' @return A GRanges object with introns removed.
#'
#' @examples
#' # Create a simple GRanges object:
#' gr <- IRanges::IRanges(
#'   start = c(13, 40, 100),
#'   end = c(20, 53, 110))
#' # Downshift it and see the introns are removed:
#' down_shift(gr)
#'
#' @export
down_shift <- function(transcript) {

  # Check if we got a GRanges object
  if (!class(transcript) %in% c("GRanges", "IRanges")) {
    stop("transcript argument must be an object of type GRanges")
  }

  for (i in 1:length(transcript)) {
    if (i == 1) {
      transcript[i] <- IRanges::shift(
        transcript[i],
        shift = 1 - start(transcript[i])
      )
    } else {
      shift_down_by <- - (start(transcript[i]) - end(transcript[i - 1])) + 1
      transcript[i] <- IRanges::shift(
        transcript[i],
        shift = shift_down_by
      )
    }
  }

  transcript
}

# -----------------------------------------------------------------------------
# Functions that validate parameters passed to functions in chimeraviz

.is_fusion_valid <- function(argument_checker, fusion) {
  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    ArgumentCheck::addError(
      msg = "'fusion' argument must be an object of type Fusion",
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is_edb_valid <- function(argument_checker, edb, fusion) {
  if (!is.null(edb)) {
    # If we got an edb object, check its validity
    if (class(edb) != "EnsDb") {
      ArgumentCheck::addError(
        msg = "'edb' argument must be an object of type EnsDb",
        argcheck = argument_checker
      )
    }
  } else {
    # If edb is not given then the fusion should have transcripts for both genes
    if (isEmpty(fusion@gene_upstream@transcripts)) {
      ArgumentCheck::addError(
        msg = paste0("There are no transcipts for gene A. Please provide an ",
                     "EnsDb object to the edb parameter, or see",
                     "get_transcripts_ensembl_db()."),
        argcheck = argument_checker
      )
    } else if (isEmpty(fusion@gene_downstream@transcripts)) {
      ArgumentCheck::addError(
        msg = paste0("There are no transcipts for gene B. Please provide an ",
                     "EnsDb object to the edb parameter, or see",
                     "get_transcripts_ensembl_db()."),
        argcheck = argument_checker
      )
    }
  }
  argument_checker
}

.is_bamfile_valid <- function(argument_checker, bamfile) {
  # Check that the argument is given
  if (is.null(bamfile) || bamfile == "") {
    ArgumentCheck::addError(
      msg = "'bamfile' must be the path to a .BAM file.",
      argcheck = argument_checker
    )
  }
  # Check that the file exists
  if (!file.exists(bamfile)) {
    ArgumentCheck::addError(
      msg = "The given 'bamfile' does not exist.",
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is_bedfile_valid <- function(argument_checker, bedfile) {
  # Check that the argument is given
  if (is.null(bedfile) || bedfile == "") {
    ArgumentCheck::addError(
      msg = "'bedfile' must be the path to a .BED file.",
      argcheck = argument_checker
    )
  }
  # Check that the file exists
  if (!file.exists(bedfile)) {
    ArgumentCheck::addError(
      msg = "The given 'bedfile' does not exist.",
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is_bedgraphfile_valid <- function(argument_checker, bedgraphfile) {
  # Check that the argument is given
  if (is.null(bedgraphfile) || bedgraphfile == "") {
    ArgumentCheck::addError(
      msg = "'bedgraphfile' must be the path to a .bedGraph file.",
      argcheck = argument_checker
    )
  }
  # Check that the file exists
  if (!file.exists(bedgraphfile)) {
    ArgumentCheck::addError(
      msg = "The given 'bedgraphfile' does not exist.",
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is_bamfile_bedgraphfile_valid <- function(
  argument_checker,
  bamfile,
  bedgraphfile) {
  # Either bamfile or bedgraphfile can be given, not both
  bamfile_given <- !is.null(bamfile)
  bedgraphfile_given <- !is.null(bedgraphfile)
  if (bamfile_given && bedgraphfile_given) {
    ArgumentCheck::addError(
      msg = "Either 'bamfile' or 'bedgraphfile' must be given, not both.",
      argcheck = argument_checker
    )
  } else if (!bamfile_given && !bedgraphfile_given) {
    # It is OK to not give any of them
  } else if (bamfile_given) {
    argument_checker <- .is_bamfile_valid(argument_checker, bamfile)
  } else {
    argument_checker <- .is_bedgraphfile_valid(argument_checker, bedgraphfile)
  }
  argument_checker
}

.is_which_transcripts_valid <- function(
  argument_checker,
  which_transcripts,
  fusion
) {
  if (class(which_transcripts) != "character") {
    ArgumentCheck::addError(
      msg = paste0("'which_transcripts' must be a character (or a character ",
                   "vector) holding the desired transcript category or the ",
                   "names of specific transcripts."),
      argcheck = argument_checker
    )
  }
  # Is which_transcripts valid?
  transcript_categories <- c(
    "exonBoundary",
    "withinExon",
    "withinIntron",
    "intergenic"
  )
  # Check if the transcript(s) given actually exist, either in
  # transcript_categories, in geneA, or in geneB
  for (i in 1:length(which_transcripts)) {
    if (
      !which_transcripts[[i]] %in% transcript_categories &&
      !which_transcripts[[i]] %in% names(fusion@gene_upstream@transcripts) &&
      !which_transcripts[[i]] %in% names(fusion@gene_downstream@transcripts)
    ) {
      ArgumentCheck::addError(
        msg = paste0("No transcript with name ", which_transcripts[[i]],
                     " was found."),
        argcheck = argument_checker
      )
    }
  }
  argument_checker
}

.is_ylim_valid <- function(argument_checker, ylim) {
  if (class(ylim) != "numeric" || length(ylim) != 2) {
    ArgumentCheck::addError(
      msg = "'ylim' must be a numeric vector of length 2",
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is_parameter_boolean <- function(
  argument_checker,
  parameter,
  parameter_name
) {
  if (class(parameter) != "logical") {
    ArgumentCheck::addError(
      msg = paste0("'", parameter_name, "'", " must be a boolean."),
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is_character_parameter_valid <- function(
  argument_checker,
  parameter,
  parameter_name
) {
  if (class(parameter) != "character") {
    ArgumentCheck::addError(
      msg = paste0("'", parameter_name, "'", " must be a character vector of ",
                   "length 1 (meaning that it's just a single string)."),
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is_nucleotide_amount_valid <- function(
  argument_checker,
  nucleotide_amount,
  fusion
) {

  fusion_junction_sequence_length <-
    length(fusion@gene_upstream@junction_sequence) +
    length(fusion@gene_downstream@junction_sequence)

  if (class(nucleotide_amount) != "numeric" ||
      nucleotide_amount <= 0 ||
      nucleotide_amount > fusion_junction_sequence_length) {
    ArgumentCheck::addError(
      msg = paste0(
        "'nucleotide_amount' must be a numeric bigger than or equal ",
        "to 0 and less than or equal to the fusion junction ",
        "sequence length."
      ),
      argcheck = argument_checker
    )
  }

  argument_checker
}

# End of functions that validate parameters passed to functions in chimeraviz
# -----------------------------------------------------------------------------

.get_transcripts_if_not_there <- function(fusion, edb) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()
  # Check parameters
  argument_checker <- .is_fusion_valid(argument_checker, fusion)
  argument_checker <- .is_edb_valid(argument_checker, edb, fusion)
  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)

  if (
    isEmpty(fusion@gene_upstream@transcripts) ||
    isEmpty(fusion@gene_downstream@transcripts)
  ) {
    message("Fetching transcripts for gene partners..")
    fusion <- get_transcripts_ensembl_db(fusion, edb)
    message("..transcripts fetched.")
  }
  fusion
}
