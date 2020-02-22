context("Test plot_fusion_transcript_with_protein_domain()")

# Load data and example fusion event
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 1)
fusion <- get_fusion_by_id(fusions, 5267)
# Select transcripts
gene_upstream_transcript <- "ENST00000434290"
gene_downstream_transcript <- "ENST00000370031"
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
# bedfile with protein domains for the transcripts in this example
bedfile <- system.file(
  "extdata",
  "protein_domains_5267.bed",
  package = "chimeraviz")
# Temporary file to store the plot
png_filename <- tempfile(
  pattern = "fusionPlot",
  fileext = ".png",
  tmpdir = tempdir())

test_that("plot_fusion_transcript_with_protein_domain produces a png file when bamfile is not given", {
  # Open device
  png(png_filename, width = 500, height = 500)
  # Plot and expect certain message output
  expect_message(
    plot_fusion_transcript_with_protein_domain(
      fusion = fusion,
      edb = edb,
      bamfile = NULL,
      bedfile = bedfile,
      gene_upstream_transcript = gene_upstream_transcript,
      gene_downstream_transcript = gene_downstream_transcript,
      plot_downstream_protein_domains_if_fusion_is_out_of_frame = TRUE
    ),
    paste0(
      "No protein domains are retained in the downstream gene, because the ",
      "fusion is not in frame. But since ",
      "plot_downstream_protein_domains_if_fusion_is_out_of_frame = TRUE, ",
      "we'll plot the protein domains in the downstream gene anyway.."
    )
  )
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})
