context("Test plotting fusion reads")

# Load data
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 1)
# Find the specific fusion we have aligned reads for
fusion <- getFusionById(fusions, 5267)
bamfile <- system.file(
  "extdata",
  "5267readsAligned.bam",
  package="chimeraviz")
# Temporary file to store the plot
pngFilename <- tempfile(
  pattern = "fusionPlot",
  fileext = ".png",
  tmpdir = tempdir())

test_that("plotFusionReads produces a png file", {
  # Expect error when we haven't added the fusion reads bam file to the fusion
  # object yet
  expect_error(plotFusionReads(fusion, pngFilename))

  # Add the bam file of aligned fusion reads to the fusion object
  fusion <- addFusionReadsAlignment(fusion, bamfile)
  # Calculate image size based on supporting reads and lenght of junction
  # sequence.
  imageWidth <- (nchar(fusion@geneA@junctionSequence) + nchar(fusion@geneB@junctionSequence)) * 10
  imageHeight <- (fusion@splitReadsCount+fusion@spanningReadsCount) * 20
  # Open device
  png(pngFilename, width = imageWidth, height = imageHeight)
  # Now we can plot
  plotFusionReads(fusion)
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(pngFilename)$size > 0, TRUE)
})
