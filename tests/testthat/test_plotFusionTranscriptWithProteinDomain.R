context("Test plotFusionTranscriptWithProteinDomain()")

# Load data and example fusion event
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 1)
fusion <- getFusionById(fusions, 5267)
# Select transcripts
geneAtranscript <- "ENST00000434290"
geneBtranscript <- "ENST00000370031"
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
# bedfile with protein domains for the transcripts in this example
bedfile <- system.file(
  "extdata",
  "protein_domains_5267.bed",
  package="chimeraviz")
# Temporary file to store the plot
pngFilename <- tempfile(
  pattern = "fusionPlot",
  fileext = ".png",
  tmpdir = tempdir())

test_that("plotFusionTranscriptWithProteinDomain produces a png file", {
  # Open device
  png(pngFilename, width = 500, height = 500)
  # Plot and expect certain message output
  expect_message(
    plotFusionTranscriptWithProteinDomain(
      fusion = fusion,
      edb = edb,
      bamfile = bamfile5267,
      bedfile = bedfile,
      geneAtranscript = geneAtranscript,
      geneBtranscript = geneBtranscript,
      plotDownstreamProteinDomainsIfFusionIsOutOfFrame = TRUE),
    "No protein domains are retained in the downstream gene, because the fusion is not in frame. But since plotDownstreamProteinDomainsIfFusionIsOutOfFrame = TRUE, we'll plot the protein domains in the downstream gene anyway..")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(pngFilename)$size > 0, TRUE)
})
