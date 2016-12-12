context("Test plotFusion")

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

test_that("plotFusion produces a png file", {
  # Temporary file to store the plot
  pngFilename <- tempfile(
    pattern = "fusionPlot",
    fileext = ".png",
    tmpdir = tempdir())
  # Open device
  png(pngFilename, width = 1000, height = 750)
  # Plot and expect certain message output
  expect_message(
    plotFusion(
      fusion = fusion,
      edb = edb,
      bamfile = bamfile5267),
    "As HENMT1 is on the minus strand, the plot for this gene will be reversed")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(pngFilename)$size > 0, TRUE)
})

test_that("plotFusion works when adding transcripts in advance", {
  # Get transcripts
  fusion <- getTranscriptsEnsembldb(fusion, edb)
  # Temporary file to store the plot
  pngFilename <- tempfile(
    pattern = "fusionPlot",
    fileext = ".png",
    tmpdir = tempdir())
  # Open device
  png(pngFilename, width = 1000, height = 750)
  # Plot and expect certain message output
  expect_message(
    plotFusion(
      fusion = fusion,
      bamfile = bamfile5267),
    "As HENMT1 is on the minus strand, the plot for this gene will be reversed")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(pngFilename)$size > 0, TRUE)
})

test_that("selectTranscript picks the correct transcript type", {
  # Expect error when we haven't added transcript to the genepartner yet
  expect_error(selectTranscript(fusion@geneA))
  # Get transcripts
  fusion <- getTranscriptsEnsembldb(fusion, edb)

  transcriptsA <- selectTranscript(fusion@geneA)
  # This GRangesList object should have 4 transcripts
  expect_equal(length(transcriptsA), 4)

  transcriptsB <- selectTranscript(fusion@geneB)
  # This GRangesList object should have 4 transcripts
  expect_equal(length(transcriptsB), 4)
})
