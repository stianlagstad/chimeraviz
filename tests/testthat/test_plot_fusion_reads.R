context("Test plot_fusion_reads()")

# Load data
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 1)
# Find the specific fusion we have aligned reads for
fusion <- get_fusion_by_id(fusions, 5267)
bamfile <- system.file(
  "extdata",
  "5267readsAligned.bam",
  package = "chimeraviz")
# Temporary file to store the plot
png_filename <- tempfile(
  pattern = "fusionPlot",
  fileext = ".png",
  tmpdir = tempdir())

test_that("plot_fusion_reads produces a png file", {
  # Expect error when we haven't added the fusion reads bam file to the fusion
  # object yet
  expect_error(plot_fusion_reads(fusion, png_filename))

  # Add the bam file of aligned fusion reads to the fusion object
  fusion <- add_fusion_reads_alignment(fusion, bamfile)
  # Calculate image size based on supporting reads and lenght of junction
  # sequence.
  image_width <- (nchar(fusion@gene_upstream@junction_sequence) +
                    nchar(fusion@gene_downstream@junction_sequence)) * 10
  image_height <- (fusion@split_reads_count + fusion@spanning_reads_count) * 20
  # Open device
  png(png_filename, width = image_width, height = image_height)
  # Now we can plot
  plot_fusion_reads(fusion)
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})
