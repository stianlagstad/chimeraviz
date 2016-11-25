context("Test utility functions")

# Import defuse results
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 1)

# Find the fusion we have get_reads.pl output for
fusion <- getFusionById(fusions, 5267)

# The read ids we're interested in
reads <- c(
  "13422259", "19375605", "29755061",
  "31632876", "32141428", "33857245")
# fastq files that has the supporting reads
fastq1 <- system.file("extdata", "reads.1.fq", package="chimeraviz")
fastq2 <- system.file("extdata", "reads.2.fq", package="chimeraviz")
# Temporary files to put extracted reads into
fastqFileOut1 <- tempfile(pattern = "fq1", tmpdir = tempdir())
fastqFileOut2 <- tempfile(pattern = "fq2", tmpdir = tempdir())

# edb object containing transcript data for the fusion transcripts from the
# deFuse example data with cluster_id=5267 and cluster_id=11759
edbSqliteFile <- system.file(
  "extdata",
  "Homo_sapiens.GRCh37.74.sqlite",
  package="chimeraviz")
edb <- ensembldb::EnsDb(edbSqliteFile)

# Get all exons for all transcripts in the genes in the fusion transcript
allTranscripts <- ensembldb::exonsBy(
  edb,
  filter = list(
    ensembldb::GeneidFilter(
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

test_that("fetchReadsFromFastq works as it should", {

  skip_on_os(os = "windows")

  fetchReadsFromFastq(reads,
                      fastq1,
                      fastq2,
                      fastqFileOut1,
                      fastqFileOut2)

  # fastqFileOut1 and fastqFileOut2 should now contain the reads supporting
  # fusion with cluster_id 5267 To test this, read the files.

  newFastqContent1 <- readr::read_lines(fastqFileOut1)
  newFastqContent2 <- readr::read_lines(fastqFileOut2)

  # We know for a fact that there are 6 reads supporting the fusion. Since each
  # read will have four lines in the file, check if the length matches
  expect_equal(length(newFastqContent1), 6 * 4)
  expect_equal(length(newFastqContent2), 6 * 4)
})

test_that("fetchReadsFromFastq fails when it's supposed to", {

  skip_on_os(os = "windows")

  expect_error(fetchReadsFromFastq(c(1,2,3), "fq1", "fq2", "fq1", "fq2"))
  expect_error(fetchReadsFromFastq(reads,
                                   fastq1,
                                   "this-is-not-an-existing-file",
                                   fastqFileOut1,
                                   fastqFileOut2))

  expect_warning(fetchReadsFromFastq(c("123123123",
                                       "213987987321",
                                       "123768"),
                                     fastq1,
                                     fastq2,
                                     fastqFileOut1,
                                     fastqFileOut2))
})

test_that("writeFusionReference writes the fusion sequence to a fasta file", {
  # Create temporary file to hold the fusion sequence
  fastaFileOut <- tempfile(pattern = "fq1", tmpdir = tempdir())
  # Write fusion sequence to file
  writeFusionReference(fusion, fastaFileOut)
  # Read the file
  fastaFileContent <- Biostrings::readDNAStringSet(fastaFileOut)

  # The fusion sequence wasn't written directly to a file, but put into a
  # DNAStringSet object first
  fusionSequence <- Biostrings::DNAStringSet(
    x = c(fusion@geneA@junctionSequence,
          fusion@geneB@junctionSequence))
  names(fusionSequence) <- "chrNA"

  # Check if equal
  expect_equal(fusionSequence, fastaFileContent)
})

test_that("writeFusionReference fails when given invalid input", {
  expect_error(writeFusionReference("not-a-fusion-object", fastaFileOut))
})

test_that("addFusionReadsAlignment works as expected", {
  # We will test this by checking the value of
  # fusion@fusionReadsAlignment@genome before and after adding the bamfile.

  # Before:
  expect_equal(is.na(fusion@fusionReadsAlignment@genome), TRUE)

  # Find alignment data
  bamfile5267 <- system.file(
    "extdata",
    "5267readsAligned.bam",
    package="chimeraviz")
  # Add the bam file of aligned fusion reads to the fusion object
  fusion <- addFusionReadsAlignment(fusion, bamfile5267)

  # After:
  expect_equal(fusion@fusionReadsAlignment@genome, fusion@genomeVersion)
})

test_that("getTranscriptsEnsembldb adds transcript data to fusion object", {
  # Add transcripts data to fusion object
  fusion <- getTranscriptsEnsembldb(fusion, edb)

  # See that we got the right amount of transcripts
  expect_equal(length(fusion@geneA@transcripts), 12)
  expect_equal(length(fusion@geneB@transcripts), 6)
})
test_that("getTranscriptsEnsembldb sets the strand if missing for a gene", {
  # This test exists because the JAFFA results does not contain strand
  # information

  fusion@geneA@strand = "*"
  fusion@geneB@strand = "*"
  fusion <- getTranscriptsEnsembldb(fusion, edb)
  expect_equal(fusion@geneA@strand, "+")
  expect_equal(fusion@geneB@strand, "-")
})

test_that("decideTranscriptCategory returns the correct categories", {
  expect_equal(decideTranscriptCategory(allTranscripts[[1]], fusion), "exonBoundary")
  expect_equal(decideTranscriptCategory(allTranscripts[[3]], fusion), "intergenic")
})

test_that("splitOnUtrAndAddFeature splits ranges as it should", {
  # Check length on the GRanges object before splitting
  expect_equal(length(allTranscripts[[1]]), 9)
  # And check it after
  expect_equal(length(splitOnUtrAndAddFeature(allTranscripts[[1]])), 11)

  # Check length on the GRanges object before splitting
  expect_equal(length(allTranscripts[[2]]), 8)
  # And check it after
  expect_equal(length(splitOnUtrAndAddFeature(allTranscripts[[2]])), 10)
})

test_that("getEnsemblIds gets the correct Ensembl gene ids for a fusion object", {
  # Store the correct Ensembl ids temporarily
  correctEnsemblIdA <- fusion@geneA@ensemblId
  correctEnsemblIdB <- fusion@geneB@ensemblId

  # Reset the fusion objects ensembl ids
  fusion@geneA@ensemblId <- ""
  fusion@geneB@ensemblId <- ""

  # Get the ensembl ids
  fusion <- getEnsemblIds(fusion)

  # Check that they are the same as the one we had originally
  expect_equal(fusion@geneA@ensemblId, correctEnsemblIdA)
  expect_equal(fusion@geneB@ensemblId, correctEnsemblIdB)

})

test_that("getEnsemblIds gives error when it can't find ensembl id", {
  # Store the correct gene names temporarily
  nameA <- fusion@geneA@name
  nameB <- fusion@geneB@name

  # Set the geneNames[1] to a non-existant gene name
  fusion@geneA@name <- "non-existant"
  fusion@geneB@name <- "non-existant"

  # Check that we get a warning
  expect_error(getEnsemblIds(fusion))

  # Restore original gene names
  fusion@geneA@name <- nameA
  fusion@geneB@name <- nameB
})

test_that("fusionToDataFrame works as expected", {
  dfFusion <- fusionToDataFrame(fusion)
  expect_is(dfFusion, "data.frame")
  expect_equal(fusion@geneA@name, as.character(dfFusion$geneA))

  # Expect error when given something other than a Fusion object
  expect_error(fusionToDataFrame(42))
})
