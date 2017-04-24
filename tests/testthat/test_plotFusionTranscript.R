context("Test plotFusionTranscript()")

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
# Temporary file to store the plot
pngFilename <- tempfile(
  pattern = "fusionPlot",
  fileext = ".png",
  tmpdir = tempdir())

test_that("plotFusionTranscripts produces a png file", {
  # Open device
  png(pngFilename, width = 500, height = 500)
  # Plot and expect certain message output
  expect_message(
    plotFusionTranscript(
      fusion = fusion,
      edb = edb),
    "Fetching transcripts for gene partners..")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(pngFilename)$size > 0, TRUE)
})
