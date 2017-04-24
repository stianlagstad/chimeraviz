context("Test selectTranscript")

# Load data and example fusion event
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 1)
fusion <- getFusionById(fusions, 5267)
# Load edb
edbSqliteFile <- system.file(
  "extdata",
  "Homo_sapiens.GRCh37.74.sqlite",
  package="chimeraviz")
edb <- ensembldb::EnsDb(edbSqliteFile)
# bamfile with reads in the regions of this fusion event
bamfile5267 <- system.file(
  "extdata",
  "fusion5267and11759reads.bam",
  package="chimeraviz")
# Get transcripts
fusion <- getTranscriptsEnsembldb(fusion, edb)

testthat::test_that("selectTranscript returns expected results with whichTranscripts=\"exonBoundary\"", {
  transcripts = selectTranscript(fusion@geneA, "exonBoundary")
  testthat::expect_equal(length(transcripts), 4)
})

testthat::test_that("selectTranscript gives exonBoundary transcripts when we ask for withinExon transcripts, because there are no witinExon transcripts", {
  testthat::expect_message(selectTranscript(fusion@geneA, "withinExon"), "..found transcripts of type exonBoundary")
})

testthat::test_that("selectTranscript gives error with whichTranscripts=\"4242\" when there are no such transcripts", {
  testthat::expect_error(selectTranscript(fusion@geneA, "4242"))
})

testthat::test_that("selectTranscript returns the expected specific transcripts", {
  transcripts <- selectTranscript(fusion@geneA, c("ENST00000373831", "ENST00000373832"))
  testthat::expect_equal(length(transcripts), 2)
})
