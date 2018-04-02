context("Test select_transcript")

# Load data and example fusion event
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 1)
fusion <- get_fusion_by_id(fusions, 5267)
# Load edb
edb_sqlite_file <- system.file(
  "extdata",
  "Homo_sapiens.GRCh37.74.sqlite",
  package = "chimeraviz")
edb <- ensembldb::EnsDb(edb_sqlite_file)
# bamfile with reads in the regions of this fusion event
bamfile5267 <- system.file(
  "extdata",
  "fusion5267and11759reads.bam",
  package = "chimeraviz")
# Get transcripts
fusion <- get_transcripts_ensembl_db(fusion, edb)

testthat::test_that(
  paste0(
    "select_transcript returns expected results with which_transcripts=",
    "\"exonBoundary\""
  ), {
    transcripts <- select_transcript(fusion@gene_upstream, "exonBoundary")
    testthat::expect_equal(length(transcripts), 4)
  }
)

testthat::test_that(
  paste0(
    "select_transcript gives exonBoundary transcripts when we ask for ",
    "withinExon transcripts, because there are no witinExon transcripts"
  )
  , {
  testthat::expect_message(
    select_transcript(fusion@gene_upstream, "withinExon"),
    "..found transcripts of type exonBoundary"
  )
})

testthat::test_that(
  paste0(
    "select_transcript gives error with which_transcripts=\"4242\" when there",
    " are no such transcripts"
  ), {
  testthat::expect_error(select_transcript(fusion@gene_upstream, "4242"))
})

testthat::test_that(
  "select_transcript returns the expected specific transcripts", {
  transcripts <- select_transcript(
    fusion@gene_upstream,
    c("ENST00000373831", "ENST00000373832")
  )
  testthat::expect_equal(length(transcripts), 2)
})
