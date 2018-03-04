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
#' @param fastqFileIn1 First fastq file to search in.
#' @param fastqFileIn2 Second fastq file to seach in.
#'
#' @param fastqFileOut1 First fastq file with results.
#' @param fastqFileOut2 Second fastq file with results.
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
#' fetchReadsFromFastq(reads, fastq1, fastq2,
#'     fastqFileOut1, fastqFileOut2)
#' # We now have the reads supporting fusion 5267 in the two files.
#' }
#'
#' @export
fetchReadsFromFastq <- function(reads,
                                fastqFileIn1,
                                fastqFileIn2,
                                fastqFileOut1,
                                fastqFileOut2) {

  # Since this function use the bash function egrep to extract reads, give a
  # warning if we're running on windows
  if (Sys.info()['sysname'] == "Windows") {
    stop(paste("This function uses the bash function egrep. It looks like",
               "you're running on windows, so this function will terminate."))
  }

  if (is.vector(reads, mode = "character") == FALSE) {
    stop("reads should be a character vector of read ids")
  }

  if (file.exists(fastqFileIn1) == FALSE || file.exists(fastqFileIn2) == FALSE) {
    stop("Invalid fastq input files")
  }

  # The command below will extract the sequences for ids 11, 22, and 33 from the
  # fastq file reads.fastq

  # $ egrep -A 3 '@11|22|33' reads.fastq | sed '/^--$/d'

  # By default, the -A parameter inserts a "--" between matches, so the sed part
  # above removes that

  # First build the query with the read ids
  query <- paste("@",
                 paste(reads, collapse = '|'),
                 sep = "")

  # Then create commands for each fastq file
  command1 <- paste("egrep -A 3 '",
                    query,
                    "' ",
                    shQuote(fastqFileIn1),
                    " | sed '/^--$/d'",
                    sep = "")
  command2 <- paste("egrep -A 3 '",
                    query,
                    "' ",
                    shQuote(fastqFileIn2),
                    " | sed '/^--$/d'",
                    sep = "")

  # We now have two commands:
  # egrep -A 1 '@11|22|33' reads.1.fastq | sed '/^--$/d'
  # egrep -A 1 '@11|22|33' reads.2.fastq | sed '/^--$/d'

  # Run the command and capture output
  supportingReadsFq1 <- system(command1, intern = TRUE)
  supportingReadsFq2 <- system(command2, intern = TRUE)

  # Write these to each their own files
  write(supportingReadsFq1, file = fastqFileOut1)
  write(supportingReadsFq1, file = fastqFileOut2)

  # Check file sizes
  if (file.info(fastqFileOut1)$size <= 10) {
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
#' fusions <- importDefuse(defuse833keFiltered, "hg19", 1)
#' # Get a specific fusion
#' fusion <- getFusionById(fusions, 5267)
#' # Create temporary file to hold the fusion sequence
#' fastaFileOut <- tempfile(pattern = "fusionSequence", tmpdir = tempdir())
#' # Write fusion sequence to file
#' writeFusionReference(fusion, fastaFileOut)
#'
#' @export
writeFusionReference <- function(fusion, filename) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  # First put the fusion junction sequence in a DNAStringSet object
  fusionSequence <- Biostrings::DNAStringSet(
    x = c(fusion@geneA@junctionSequence,
          fusion@geneB@junctionSequence))

  # Give an error if the length of the fusionSequence is 0:
  if (nchar(fusionSequence) == 0) {
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
  names(fusionSequence) <- "chrNA"
  # Write to file
  Biostrings::writeXStringSet(fusionSequence, filename)
}

#' Get ensembl ids for a fusion object
#'
#' This function will get the ensembl ids from the org.Hs.eg.db package given
#' the gene names of the fusion event.
#'
#' @param fusion The Fusion object we want to get ensembl ids for.
#'
#' @return The Fusion object with Ensembl ids set.
#'
#' @import org.Hs.eg.db
#'
#' @examples
#' # Import the filtered defuse results
#' defuse833keFiltered <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833keFiltered, "hg19", 1)
#' # Get a specific fusion
#' fusion <- getFusionById(fusions, 5267)
#' # See the ensembl ids:
#' partnerGeneEnsemblId(upstreamPartnerGene(fusion))
#' # [1] "ENSG00000180198"
#' partnerGeneEnsemblId(downstreamPartnerGene(fusion))
#' # [1] "ENSG00000162639"
#' # Reset the fusion objects ensembl ids
#' partnerGeneEnsemblId(upstreamPartnerGene(fusion)) <- ""
#' partnerGeneEnsemblId(downstreamPartnerGene(fusion))  <- ""
#' # Get the ensembl ids
#' fusion <- getEnsemblIds(fusion)
#' # See that we now have the same ensembl ids again:
#' partnerGeneEnsemblId(upstreamPartnerGene(fusion))
#' # [1] "ENSG00000180198"
#' partnerGeneEnsemblId(downstreamPartnerGene(fusion))
#' # [1] "ENSG00000162639"
#'
#' @export
getEnsemblIds <- function(fusion) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  result <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = c(fusion@geneA@name, fusion@geneB@name),
    keytype = "ALIAS",
    columns = c("ALIAS", "ENSEMBL"))

  geneAresults <- result[which(result$ALIAS == fusion@geneA@name), ]
  geneBresults <- result[which(result$ALIAS == fusion@geneB@name), ]

  # Stop execution if no results
  if (any(is.na(geneAresults$ENSEMBL))) {
    stop(paste("Could not find Ensembl id for ", fusion@geneA@name, ". ",
               "If you know the id, add it manually with ",
               "fusion@geneA@ensemblId <- \"ensemblId\"",
               sep = ""))
  }
  if (any(is.na(geneBresults$ENSEMBL))) {
    stop(paste("Could not find Ensembl id for ", fusion@geneB@name, ". ",
               "If you know the id, add it manually with ",
               "fusion@geneB@ensemblId <- \"ensemblId\"",
               sep = ""))
  }

  # Store ensembl ids in fusion object
  fusion@geneA@ensemblId <- geneAresults[, 2]
  fusion@geneB@ensemblId <- geneBresults[, 2]

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
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- getFusionById(fusions, 5267)
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
#'         partnerGeneEnsemblId(upstreamPartnerGene(fusion)),
#'         partnerGeneEnsemblId(downstreamPartnerGene(fusion))))),
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
#' gr <- splitOnUtrAndAddFeature(gr)
#' # Check the length again
#' length(gr)
#' # Should be 11 now, as the range containing the cds_strat position and the
#' # range containing the cds_stop position has been split into separate ranges
#'
#' @importFrom S4Vectors mcols
#'
#' @export
splitOnUtrAndAddFeature <- function(gr) {

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
    end(first) <- cds_start-1
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
    start(second) <- cds_end+1
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
  #     end(gr) < cds_start
  #   minus strand:
  #     start(gr) > cds_end
  if (length(gr[as.character(strand(gr)) == "+" & end(gr) < cds_start |
                as.character(strand(gr)) == "-" & start(gr) > cds_end])) {
    S4Vectors::mcols(gr[as.character(strand(gr)) == "+" & end(gr) < cds_start |
               as.character(strand(gr)) == "-" & start(gr) > cds_end])$feature <- "5utr"
  }
  # 3utr:
  #   plus strand:
  #     start(gr) > cds_end
  #   minus strand:
  #     end(gr) < cds_start
  if (length(gr[as.character(strand(gr)) == "+" & start(gr) > cds_end |
                as.character(strand(gr)) == "-" & end(gr) < cds_start])) {
    S4Vectors::mcols(gr[as.character(strand(gr)) == "+" & start(gr) > cds_end |
               as.character(strand(gr)) == "-" & end(gr) < cds_start])$feature <- "3utr"
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
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- getFusionById(fusions, 5267)
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
#'         partnerGeneEnsemblId(upstreamPartnerGene(fusion)),
#'         partnerGeneEnsemblId(downstreamPartnerGene(fusion))))),
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
#' decideTranscriptCategory(gr, fusion)
#' # "exonBoundary"
#' # Check another case
#' gr <- allTranscripts[[3]]
#' decideTranscriptCategory(gr, fusion)
#' # "withinIntron"
#'
#' @importFrom S4Vectors mcols
#'
#' @export
decideTranscriptCategory <- function(gr, fusion) {

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
  geneA <- S4Vectors::mcols(gr)$gene_id[[1]] == fusion@geneA@ensemblId
  breakpoint <- if(all(S4Vectors::mcols(gr)$gene_id == fusion@geneA@ensemblId)) fusion@geneA@breakpoint else fusion@geneB@breakpoint
  exonStartPositions <- GenomicRanges::start(gr)
  exonEndPositions <- GenomicRanges::end(gr)

  # Check whether or not the breakpoint is within the transcript
  if (all(breakpoint < exonStartPositions) || all(exonEndPositions < breakpoint)) {
    return("intergenic")
  }

  # Assuming only single-stranded transcripts
  if (as.character(strand(gr)[1]) == "+") {

    # To decide whether or not the breakpoint occurs at an exon boundary, we
    # have to be vary careful about both whether this is the upstream/downstrem
    # fusion partner gene
    if (geneA) {
      # GeneA
      if (any(exonEndPositions == breakpoint)) {
        # exon boundary
        return("exonBoundary")
      }
    } else {
      # GeneB
      if (any(exonStartPositions == breakpoint)) {
        # exon boundary
        return("exonBoundary")
      }
    }

  } else {

    # To decide whether or not the breakpoint occurs at an exon boundary, we
    # have to be vary careful about both whether this is the upstream/downstrem
    # fusion partner gene
    if (all(geneA)) {
      # GeneA
      if (any(exonStartPositions == breakpoint)) {
        # exon boundary
        return("exonBoundary")
      }
    } else {
      # GeneB
      if (any(exonEndPositions == breakpoint)) {
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
#' updated Fusion object, where the fusion@geneA@transcriptsX slots are set with
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
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- getFusionById(fusions, 5267)
#' # Create edb object
#' edbSqliteFile <- system.file(
#'   "extdata",
#'   "Homo_sapiens.GRCh37.74.sqlite",
#'   package="chimeraviz")
#' edb <- ensembldb::EnsDb(edbSqliteFile)
#' # Add transcripts data to fusion object
#' fusion <- getTranscriptsEnsembldb(fusion, edb)
#' # The transcripts are now accessible through fusion@geneA@transcripts and
#' # fusion@geneB@transcripts .
#'
#' @importFrom S4Vectors mcols
#' @importFrom ensembldb exonsBy
#' @importFrom AnnotationFilter GeneIdFilter
#'
#' @export
getTranscriptsEnsembldb <- function(fusion, edb) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  # Check if we got a valid edb
  if (class(edb) != "EnsDb") {
    stop("edb argument must be an object of type EnsDb")
  }

  # Fetch ensembl gene ids if we don't already have them
  if (is.na(fusion@geneA@ensemblId) || is.na(fusion@geneB@ensemblId)) {
    fusion <- getEnsemblIds(fusion)
  }

  # Get all exon information
  allTranscripts <- ensembldb::exonsBy(
    edb,
    filter = list(
      AnnotationFilter::GeneIdFilter(
        list(
          fusion@geneA@ensemblId,
          fusion@geneB@ensemblId))),
    columns = c(
      "gene_id",
      "gene_name",
      "tx_id",
      "tx_cds_seq_start",
      "tx_cds_seq_end",
      "exon_id"))

  # Fail if no transcripts were found
  if (length(allTranscripts) == 0) {
    stop(paste(
      "No transcripts available for the genes ",
      fusion@geneA@name,
      " and ",
      fusion@geneB@name,
      ".",
      sep = ""))
  }

  # Go through each transcript in the GRangesList and
  grlA <- GRangesList()
  grlB <- GRangesList()
  for(i in 1:length(allTranscripts)) {

    # Extract the GRanges object
    gr <- allTranscripts[[i]]

    # Add $transcript as a metadata column
    S4Vectors::mcols(gr)$transcript <- names(allTranscripts[i])
    # Add $symbol as a metadata column
    S4Vectors::mcols(gr)$symbol <- names(allTranscripts[i])

    # Add utr/coding S4Vectors::mcols$feature data for each region
    if(is.na(S4Vectors::mcols(gr)$tx_cds_seq_start[[1]])) {
      # The transcript is not coding. Set all to "utr"
      S4Vectors::mcols(gr)$feature <- "utr"
    }else {
      # The transcript is coding. The function below will set
      # S4Vectors::mcols(gr)$feature = "utr5"/"protein_coding"/"utr3" for each
      # region
      gr <- splitOnUtrAndAddFeature(gr)
    }

    # Create new GRangesList
    grl <- GRangesList(gr)

    # Tag each transcript with either "exonBoundary", "withinExon",
    # "withinIntron", or "intergenic", depending on where in the transcript the
    # fusion breakpoint hits
    S4Vectors::mcols(grl)$transcriptCategory <- decideTranscriptCategory(
      gr,
      fusion)

    # Set S4Vectors::mcols$transcript and name for the grl
    S4Vectors::mcols(grl)$transcript <- names(allTranscripts[i])
    names(grl) <- names(allTranscripts[i])

    # Add S4Vectors::mcols$coding to signify whether it is a coding transcript
    # at all
    S4Vectors::mcols(grl)$coding <-
      !all(is.na(S4Vectors::mcols(gr)$tx_cds_seq_start))

    # Append the GRanges object to the correct GRangesList
    if (S4Vectors::mcols(gr)$gene_id[[1]] == fusion@geneA@ensemblId) {
      grlA <- append(grlA, grl)
    } else {
      grlB <- append(grlB, grl)
    }
  }

  # Add transcripts to fusion object
  fusion@geneA@transcripts <- grlA
  fusion@geneB@transcripts <- grlB

  # Warn if no transcripts were found for one of the genes
  if (length(fusion@geneA@transcripts) == 0) {
    warning(paste(
      "No transcripts available for the upstream gene ",
      fusion@geneA@name,
      " available.",
      sep = ""))
  }
  if (length(fusion@geneA@transcripts) == 0) {
    warning(paste(
      "No transcripts available for the downstream gene ",
      fusion@geneB@name,
      " available.",
      sep = ""))
  }

  # In case the fusion object doesn't have the strands set (as is the case for
  # JAFFA), set the strands now:
  if (fusion@geneA@strand == "*" || fusion@geneB@strand == "*") {
    fusion@geneA@strand = as.character(strand(fusion@geneA@transcripts[[1]][1]))
    fusion@geneB@strand = as.character(strand(fusion@geneB@transcripts[[1]][1]))
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
#' @return An updated fusion object with fusion@fusionReadsAlignment set.
#'
#' @examples
#' # Load data
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 1)
#' # Find the specific fusion we have aligned reads for
#' fusion <- getFusionById(fusions, 5267)
#' # Get reference to the bamfile with the alignment data
#' bamfile5267 <- system.file(
#'   "extdata",
#'   "5267readsAligned.bam",
#'   package="chimeraviz")
#' # Add the bam file of aligned fusion reads to the fusion object
#' fusion <- addFusionReadsAlignment(fusion, bamfile5267)
#'
#' @export
addFusionReadsAlignment <- function(fusion, bamfile) {

  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    stop("fusion argument must be an object of type Fusion")
  }

  fusionReadsAlignment <- Gviz::AlignmentsTrack(
    bamfile,
    isPaired = TRUE,
    # Setting chromosome to chrNA because this is a fusion sequence not found in
    # any reference genome.
    chromosome = "chrNA",
    name="Fusion Reads",
    genome = fusion@genomeVersion)

  # Return new fusion object, now with the fusion read alignment
  fusion@fusionReadsAlignment = fusionReadsAlignment
  fusion
}

#' Coerce Fusion object to data.frame
#'
#' This function is used in createFusionReport() to convert Fusion objects to a
#' data.frame-format.
#'
#' @param fusion The Fusion object to coerce.
#'
#' @return A data.frame with the fusion object.
#'
#' @seealso createFusionReport
#'
#' @examples
#' # Load data
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 1)
#' # Find the fusion object to create a data frame from
#' fusion <- getFusionById(fusions, 5267)
#' # Create the data frame
#' dfFusion <- fusionToDataFrame(fusion)
#'
#' @export
fusionToDataFrame <- function(fusion) {

  # Check if we got a list of fusion objects
  if (class(fusion) != "Fusion") {
    stop("fusions argument must be a list of Fusion objects")
  }

  df <- data.frame(
    fusion@id,
    fusion@geneA@name,
    fusion@geneA@ensemblId,
    fusion@geneA@breakpoint,
    fusion@geneB@name,
    fusion@geneB@ensemblId,
    fusion@geneB@breakpoint,
    fusion@splitReadsCount,
    fusion@spanningReadsCount)
  names(df) <- c(
    "id",
    "geneA",
    "ensemblA",
    "breakpointA",
    "geneB",
    "ensemblB",
    "breakpointB",
    "Split Reads",
    "Spanning Reads")
  df
}

#' Select which transcript to use (for plotting) for a GenePartner object
#'
#' This function takes a GenePartner object and creates a transcript data.frame
#' with transcript information, including only the transcripts given by the
#' parameter whichTranscripts.
#'
#' selectTranscript() selects which transcript to create by this prioritization:
#'
#' 1. Exon boundary transcripts.
#' 2. Within exon transcripts.
#' 3. Within intron transcripts.
#' 4. Intergenic transcripts.
#'
#' @param genePartner The GenePartner object to select a transcript for.
#' @param whichTranscripts This character vector decides which transcripts are
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
#' fusions <- importDefuse(defuse833ke, "hg19", 1)
#' fusion <- getFusionById(fusions, 5267)
#' # Load edb
#' edbSqliteFile <- system.file(
#'   "extdata",
#'   "Homo_sapiens.GRCh37.74.sqlite",
#'   package="chimeraviz")
#' edb <- ensembldb::EnsDb(edbSqliteFile)
#' # Get transcripts
#' fusion <- getTranscriptsEnsembldb(fusion, edb)
#' # Select transcript
#' transcriptsA <- selectTranscript(upstreamPartnerGene(fusion))
#'
#' @export
selectTranscript <- function(
  genePartner,
  whichTranscripts = "exonBoundary") {

  # Check if we got a PartnerGene object
  if (class(genePartner) != "PartnerGene") {
    stop("genePartner argument must be an object of type PartnerGene")
  }

  # Does the PartnerGene have any transcripts?
  if (isEmpty(genePartner@transcripts)) {
    stop("genePartner has no transcripts. See getTranscriptsEnsembldb()")
  }

  # Either select from one of these four categories
  # 1. Exon boundary transcripts.
  # 2. Within exon transcripts.
  # 3. Within intron transcripts.
  # 4. Intergenic transcripts.
  #
  # or select the specific transcripts given by whichTranscripts.

  transcriptCategories <- c("exonBoundary", "withinExon", "withinIntron", "intergenic")

  if (whichTranscripts[[1]] %in% transcriptCategories) {

    message(paste0("Selecting transcripts for ", genePartner@name, ".."))

    # If the user has chosen one of the four transcript categories, then we want to check whether or not such
    # transcripts exist. If they exist, simply return them. If they don't exist, go on to try the other categories. Try
    # the wanted category first:
    if (length(genePartner@transcripts[mcols(genePartner@transcripts)$transcriptCategory == whichTranscripts[[1]] ]) > 0) {
      message(paste0("..found transcripts of type ", whichTranscripts[[1]]))
      return(genePartner@transcripts[mcols(genePartner@transcripts)$transcriptCategory == whichTranscripts[[1]] ])
    }
    # Check the remaining categories
    for(transcriptCategory in transcriptCategories[transcriptCategories != whichTranscripts[[1]]]) {
      if (length(genePartner@transcripts[mcols(genePartner@transcripts)$transcriptCategory == transcriptCategory ]) > 0) {
        message(paste0("..found transcripts of type ", transcriptCategory))
        return(genePartner@transcripts[mcols(genePartner@transcripts)$transcriptCategory == transcriptCategory ])
      }
    }
  }

  # At this point the user wants specific transcripts. Get the transcripts that
  # we have for this genePartner.

  if (length(genePartner@transcripts[names(genePartner@transcripts) %in% whichTranscripts]) > 0) {
    return(genePartner@transcripts[names(genePartner@transcripts) %in% whichTranscripts])
  }

  stop("The specific transcripts could not be found")
}

# Check that there's at least one transcript that has the fusion breakpoint
# within the transcript.
.checkThatBreakpointsAreWithinTranscripts <- function(
  fusion,
  transcriptsA,
  transcriptsB) {
  if (!any(start(transcriptsA) < fusion@geneA@breakpoint) &&
      any(fusion@geneA@breakpoint < end(transcriptsA))) {
    stop(paste0(
      "None of the transcripts given for gene A has the fusion breakpoint ",
      "within them. This plot cannot be created with the given transcripts."))
  }
  if (!any(start(transcriptsB) < fusion@geneB@breakpoint) &&
      any(fusion@geneB@breakpoint < end(transcriptsB))) {
    stop(paste0(
      "None of the transcripts given for gene B has the fusion breakpoint ",
      "within them. This plot cannot be created with the given transcripts."))
  }
}

# Check that the transcripts have a breakpoint exon
.checkThatTranscriptsHaveBreakpointExons <- function(
  fusion,
  transcriptA,
  transcriptB) {

  if (fusion@geneA@strand == "+") {
    if (length(transcriptA[end(ranges(transcriptA)) == fusion@geneA@breakpoint]) == 0) {
      stop(paste(
        "The transcript for gene A doesn't have an exon boundary matching the",
        "fusion breakpoint. The plot cannot be created."))
    }
  } else {
    if (length(transcriptA[start(ranges(transcriptA)) == fusion@geneA@breakpoint]) == 0) {
      stop(paste(
        "The transcript for gene A doesn't have an exon boundary matching the",
        "fusion breakpoint. The plot cannot be created."))
    }
  }
  if (fusion@geneB@strand == "+") {
    if (length(transcriptB[start(ranges(transcriptB)) == fusion@geneB@breakpoint]) == 0) {
      stop(paste(
        "The transcript for gene B doesn't have an exon boundary matching the",
        "fusion breakpoint. The plot cannot be created."))
    }
  } else {
    if (length(transcriptB[end(ranges(transcriptB)) == fusion@geneB@breakpoint]) == 0) {
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
#' downShift(gr)
#'
#' @export
downShift <- function(transcript) {

  # Check if we got a GRanges object
  if (!class(transcript) %in% c("GRanges", "IRanges")) {
    stop("transcript argument must be an object of type GRanges")
  }

  for (i in 1:length(transcript)) {
    if (i == 1) {
      transcript[i] <- IRanges::shift(transcript[i], shift = 1-start(transcript[i]))
    } else {
      shiftDownBy <- -(start(transcript[i])-end(transcript[i-1]))+1
      transcript[i] <- IRanges::shift(transcript[i], shift = shiftDownBy)
    }
  }

  transcript
}

# -----------------------------------------------------------------------------
# Functions that validate parameters passed to functions in chimeraviz

.is.fusion.valid <- function(argument_checker, fusion) {
  # Check if we got a fusion object
  if (class(fusion) != "Fusion") {
    ArgumentCheck::addError(
      msg = "'fusion' argument must be an object of type Fusion",
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is.edb.valid <- function(argument_checker, edb, fusion) {
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
    if (isEmpty(fusion@geneA@transcripts)) {
      ArgumentCheck::addError(
        msg = paste0("There are no transcipts for gene A. Please provide an ",
                     "EnsDb object to the edb parameter, or see",
                     "getTranscriptsEnsembldb()."),
        argcheck = argument_checker
      )
    } else if (isEmpty(fusion@geneB@transcripts)) {
      ArgumentCheck::addError(
        msg = paste0("There are no transcipts for gene B. Please provide an ",
                     "EnsDb object to the edb parameter, or see",
                     "getTranscriptsEnsembldb()."),
        argcheck = argument_checker
      )
    }
  }
  argument_checker
}

.is.bamfile.valid <- function(argument_checker, bamfile) {
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

.is.bedfile.valid <- function(argument_checker, bedfile) {
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

.is.bedgraphfile.valid <- function(argument_checker, bedgraphfile) {
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

.is.either.bamfile.or.bedgraphfile.valid <- function(argument_checker, bamfile, bedgraphfile) {
  # Either bamfile or bedgraphfile can be given, not both
  bamfileGiven <- !is.null(bamfile)
  bedgraphfileGiven <- !is.null(bedgraphfile)
  if (bamfileGiven && bedgraphfileGiven) {
    ArgumentCheck::addError(
      msg = "Either 'bamfile' or 'bedgraphfile' must be given, not both.",
      argcheck = argument_checker
    )
  } else if (!bamfileGiven && !bedgraphfileGiven) {
    ArgumentCheck::addError(
      msg = "Either 'bamfile' or 'bedgraphfile' must be given",
      argcheck = argument_checker
    )
  } else if (bamfileGiven) {
    argument_checker <- .is.bamfile.valid(argument_checker, bamfile)
  } else {
    argument_checker <- .is.bedgraphfile.valid(argument_checker, bedgraphfile)
  }
  argument_checker
}

.is.whichTranscripts.valid <- function(argument_checker, whichTranscripts, fusion) {
  if (class(whichTranscripts) != "character") {
    ArgumentCheck::addError(
      msg = paste0("'whichTranscripts' must be a character (or a character ",
                   "vector) holding the desired transcript category or the ",
                   "names of specific transcripts."),
      argcheck = argument_checker
    )
  }
  # Is whichTranscripts valid?
  transcriptCategories <- c(
    "exonBoundary",
    "withinExon",
    "withinIntron",
    "intergenic"
  )
  # Check if the transcript(s) given actually exist, either in
  # transcriptCategories, in geneA, or in geneB
  for (i in 1:length(whichTranscripts)) {
    if (
      !whichTranscripts[[i]] %in% transcriptCategories &&
      !whichTranscripts[[i]] %in% names(fusion@geneA@transcripts) &&
      !whichTranscripts[[i]] %in% names(fusion@geneB@transcripts)
    ) {
      ArgumentCheck::addError(
        msg = paste0("No transcript with name ", whichTranscripts[[i]],
                     " was found."),
        argcheck = argument_checker
      )
    }
  }
  argument_checker
}

.is.ylim.valid <- function(argument_checker, ylim) {
  if (class(ylim) != "numeric" || length(ylim) != 2) {
    ArgumentCheck::addError(
      msg = "'ylim' must be a numeric vector of length 2",
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is.parameter.boolean <- function(argument_checker, parameter, parameterName) {
  if (class(parameter) != "logical") {
    ArgumentCheck::addError(
      msg = paste0("'", parameterName, "'", " must be a boolean."),
      argcheck = argument_checker
    )
  }
  argument_checker
}

.is.character.parameter.valid <- function(argument_checker, parameter, parameterName) {
  if (class(parameter) != "character") {
    ArgumentCheck::addError(
      msg = paste0("'", parameterName, "'", " must be a character vector of ",
                   "length 1 (meaning that it's just a single string)."),
      argcheck = argument_checker
    )
  }
  argument_checker
}

is.nucleotideAmount.valid <- function(argument_checker, nucleotideAmount, fusion) {

  fusionJunctionSequenceLength <-
    length(fusion@geneA@junctionSequence) +
    length(fusion@geneB@junctionSequence)

  if (class(nucleotideAmount) != "numeric" ||
      nucleotideAmount <= 0 ||
      nucleotideAmount > fusionJunctionSequenceLength) {
    ArgumentCheck::addError(
      msg = paste0("'nucleotideAmount' must be a numeric bigger than or equal ",
                   "to 0 and less than or equal to the fusion junction ",
                   "sequence length."),
      argcheck = argument_checker
    )
  }

  argument_checker
}

# End of functions that validate parameters passed to functions in chimeraviz
# -----------------------------------------------------------------------------

.getTranscriptsIfNotThere <- function(fusion, edb) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()
  # Check parameters
  argument_checker <- .is.fusion.valid(argument_checker, fusion)
  argument_checker <- .is.edb.valid(argument_checker, edb, fusion)
  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)

  if (isEmpty(fusion@geneA@transcripts) || isEmpty(fusion@geneB@transcripts)) {
    message("Fetching transcripts for gene partners..")
    fusion <- getTranscriptsEnsembldb(fusion, edb)
    message("..transcripts fetched.")
  }
  fusion
}
