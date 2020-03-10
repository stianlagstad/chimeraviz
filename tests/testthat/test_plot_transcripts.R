context("Test plot_transcripts()")

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
# Temporary file to store the plot
png_filename <- tempfile(
  pattern = "fusionPlot",
  fileext = ".png",
  tmpdir = tempdir())

test_that("plot_transcripts produces a png file", {
  # Open device
  png(png_filename, width = 500, height = 500)
  # Plot and expect certain message output
  expect_message(
    plot_transcripts(
      fusion = fusion,
      edb = edb),
    "Fetching transcripts for gene partners..")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})

test_that("plot_transcripts with a bamfile produces a png file", {
  bamfile5267 <- system.file(
    "extdata",
    "fusion5267and11759reads.bam",
    package = "chimeraviz")
  # Open device
  png(png_filename, width = 500, height = 500)
  # Plot and expect certain message output
  expect_message(
    plot_transcripts(
      fusion = fusion,
      bamfile = bamfile5267,
      edb = edb),
    "Fetching transcripts for gene partners..")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})

test_that("plot_transcripts with a bedgraphfile produces a png file", {
  bedgraphfile <- system.file(
    "extdata",
    "fusion5267and11759reads.bedGraph",
    package = "chimeraviz"
  )
  # Open device
  png(png_filename, width = 500, height = 500)
  # Plot and expect certain message output
  expect_message(
    plot_transcripts(
      fusion = fusion,
      bedgraphfile = bedgraphfile,
      edb = edb),
    "Fetching transcripts for gene partners..")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})

test_that("plot_transcripts with reduce_transcripts produces a png file", {
  # Open device
  png(png_filename, width = 500, height = 500)
  # Plot and expect certain message output
  expect_message(
    plot_transcripts(
      fusion = fusion,
      reduce_transcripts = TRUE,
      edb = edb),
    "Fetching transcripts for gene partners..")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})
