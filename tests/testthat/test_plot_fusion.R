context("Test plot_fusion()")

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

test_that("plot_fusion produces a png file", {
  # Temporary file to store the plot
  png_filename <- tempfile(
    pattern = "fusionPlot",
    fileext = ".png",
    tmpdir = tempdir())
  # Open device
  png(png_filename, width = 1000, height = 750)
  # Plot and expect certain message output
  expect_message(
    plot_fusion(
      fusion = fusion,
      edb = edb,
      bamfile = bamfile5267),
    "As HENMT1 is on the minus strand, the plot for this gene will be reversed")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})

test_that("plot_fusion works when adding transcripts in advance", {
  # Get transcripts
  fusion <- get_transcripts_ensembl_db(fusion, edb)
  # Temporary file to store the plot
  png_filename <- tempfile(
    pattern = "fusionPlot",
    fileext = ".png",
    tmpdir = tempdir())
  # Open device
  png(png_filename, width = 1000, height = 750)
  # Plot and expect certain message output
  expect_message(
    plot_fusion(
      fusion = fusion,
      bamfile = bamfile5267),
    "As HENMT1 is on the minus strand, the plot for this gene will be reversed")
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})

test_that("select_transcript picks the correct transcript type", {
  # Expect error when we haven't added transcript to the genepartner yet
  expect_error(select_transcript(fusion@gene_upstream))
  # Get transcripts
  fusion <- get_transcripts_ensembl_db(fusion, edb)

  transcripts_upstream <- select_transcript(fusion@gene_upstream)
  # This GRangesList object should have 4 transcripts
  expect_equal(length(transcripts_upstream), 4)

  transcripts_downstream <- select_transcript(fusion@gene_downstream)
  # This GRangesList object should have 4 transcripts
  expect_equal(length(transcripts_downstream), 4)
})
