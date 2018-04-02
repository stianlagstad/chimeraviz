context("Test utility functions")

# Import defuse results
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 1)

# Find the fusion we have get_reads.pl output for
fusion <- get_fusion_by_id(fusions, 5267)

# The read ids we're interested in
reads <- c(
  "13422259", "19375605", "29755061",
  "31632876", "32141428", "33857245")
# fastq files that has the supporting reads
fastq1 <- system.file("extdata", "reads.1.fq", package = "chimeraviz")
fastq2 <- system.file("extdata", "reads.2.fq", package = "chimeraviz")
# Temporary files to put extracted reads into
fastq_file_out_1 <- tempfile(pattern = "fq1", tmpdir = tempdir())
fastq_file_out_2 <- tempfile(pattern = "fq2", tmpdir = tempdir())

# edb object containing transcript data for the fusion transcripts from the
# deFuse example data with cluster_id=5267 and cluster_id=11759
edb_sqlite_file <- system.file(
  "extdata",
  "Homo_sapiens.GRCh37.74.sqlite",
  package = "chimeraviz")
edb <- ensembldb::EnsDb(edb_sqlite_file)

# Get all exons for all transcripts in the genes in the fusion transcript
all_transcripts <- ensembldb::exonsBy(
  edb,
  filter = list(
    AnnotationFilter::GeneIdFilter(
      list(
        fusion@gene_upstream@ensembl_id,
        fusion@gene_downstream@ensembl_id))),
  columns = c(
    "gene_id",
    "gene_name",
    "tx_id",
    "tx_cds_seq_start",
    "tx_cds_seq_end",
    "exon_id"))

test_that("fetch_reads_from_fastq works as it should", {

  skip_on_os(os = "windows")

  fetch_reads_from_fastq(
    reads,
    fastq1,
    fastq2,
    fastq_file_out_1,
    fastq_file_out_2)

  # fastqFileOut1 and fastqFileOut2 should now contain the reads supporting
  # fusion with cluster_id 5267 To test this, read the files.

  new_fastq_content_1 <- readr::read_lines(fastq_file_out_1)
  new_fastq_content_2 <- readr::read_lines(fastq_file_out_2)

  # We know for a fact that there are 6 reads supporting the fusion. Since each
  # read will have four lines in the file, check if the length matches
  expect_equal(length(new_fastq_content_1), 6 * 4)
  expect_equal(length(new_fastq_content_2), 6 * 4)
})

test_that("fetch_reads_from_fastq fails when it's supposed to", {

  skip_on_os(os = "windows")

  expect_error(
    fetch_reads_from_fastq(
      c(1, 2, 3),
      "fq1",
      "fq2",
      "fq1",
      "fq2")
    )
  expect_error(
    fetch_reads_from_fastq(
      reads,
      fastq1,
      "this-is-not-an-existing-file",
      fastq_file_out_1,
      fastq_file_out_2)
    )

  expect_warning(
    fetch_reads_from_fastq(
      c("123123123", "213987987321", "123768"),
      fastq1,
      fastq2,
      fastq_file_out_1,
      fastq_file_out_2)
    )
})

test_that("write_fusion_reference writes the fusion sequence to a fasta file", {
  # Create temporary file to hold the fusion sequence
  fasta_file_out <- tempfile(pattern = "fq1", tmpdir = tempdir())
  # Write fusion sequence to file
  write_fusion_reference(fusion, fasta_file_out)
  # Read the file
  fasta_file_content <- Biostrings::readDNAStringSet(fasta_file_out)

  # The fusion sequence wasn't written directly to a file, but put into a
  # DNAStringSet object first
  fusion_sequence <- Biostrings::DNAStringSet(
    x = c(fusion@gene_upstream@junction_sequence,
          fusion@gene_downstream@junction_sequence))
  names(fusion_sequence) <- "chrNA"

  # Check if equal
  expect_equal(fusion_sequence, fasta_file_content)
})

test_that("write_fusion_reference fails when given invalid input", {
  expect_error(write_fusion_reference("not-a-fusion-object", fasta_file_out))
})

test_that(
  paste0(
    "write_fusion_reference fails when there's no fusion junction sequence ",
    "in the Fusion object"
  ), {
  # Create a Fusion object with no fusion junction sequences
  fusion_no_Junc_seq <- fusion
  fusion_no_Junc_seq@gene_upstream@junction_sequence <-
    Biostrings::DNAString("")
  fusion_no_Junc_seq@gene_downstream@junction_sequence <-
    Biostrings::DNAString("")
  # Create temporary file to hold the fusion sequence
  fasta_file_out <- tempfile(pattern = "fq1", tmpdir = tempdir())
  # Call the function
  expect_error(write_fusion_reference(fusion_no_Junc_seq, fasta_file_out))
})

test_that("add_fusion_reads_alignment works as expected", {
  # We will test this by checking the value of
  # fusion@fusion_reads_alignment@genome before and after adding the bamfile.

  # Before:
  expect_equal(is.na(fusion@fusion_reads_alignment@genome), TRUE)

  # Find alignment data
  bamfile5267 <- system.file(
    "extdata",
    "5267readsAligned.bam",
    package = "chimeraviz")
  # Add the bam file of aligned fusion reads to the fusion object
  fusion <- add_fusion_reads_alignment(fusion, bamfile5267)

  # After:
  expect_equal(fusion@fusion_reads_alignment@genome, fusion@genome_version)
})

test_that("get_transcripts_ensembl_db adds transcript data to fusion object", {
  # Add transcripts data to fusion object
  fusion <- get_transcripts_ensembl_db(fusion, edb)

  # See that we got the right amount of transcripts
  expect_equal(length(fusion@gene_upstream@transcripts), 12)
  expect_equal(length(fusion@gene_downstream@transcripts), 6)
})
test_that("get_transcripts_ensembl_db sets the strand if missing for a gene", {
  # This test exists because the JAFFA results does not contain strand
  # information

  fusion@gene_upstream@strand <- "*"
  fusion@gene_downstream@strand <- "*"
  fusion <- get_transcripts_ensembl_db(fusion, edb)
  expect_equal(fusion@gene_upstream@strand, "+")
  expect_equal(fusion@gene_downstream@strand, "-")
})

test_that("decide_transcript_category returns the correct categories", {
  expect_equal(
    decide_transcript_category(all_transcripts[[1]], fusion),
    "exonBoundary"
  )
  expect_equal(
    decide_transcript_category(all_transcripts[[3]], fusion),
    "intergenic"
  )
})

test_that("split_on_utr_and_add_feature splits ranges as it should", {
  # Check length on the GRanges object before splitting
  expect_equal(length(all_transcripts[[1]]), 9)
  # And check it after
  expect_equal(length(split_on_utr_and_add_feature(all_transcripts[[1]])), 11)

  # Check length on the GRanges object before splitting
  expect_equal(length(all_transcripts[[2]]), 8)
  # And check it after
  expect_equal(length(split_on_utr_and_add_feature(all_transcripts[[2]])), 10)
})

test_that(
  "get_ensembl_ids gets the correct Ensembl gene ids for a fusion object", {
  # Store the correct Ensembl ids temporarily
  correct_ensembl_id_upstream <- fusion@gene_upstream@ensembl_id
  correct_ensembl_id_downstream <- fusion@gene_downstream@ensembl_id

  # Reset the fusion objects ensembl ids
  fusion@gene_upstream@ensembl_id <- ""
  fusion@gene_downstream@ensembl_id <- ""

  # Get the ensembl ids
  fusion <- get_ensembl_ids(fusion)

  # Check that they are the same as the one we had originally
  expect_equal(fusion@gene_upstream@ensembl_id, correct_ensembl_id_upstream)
  expect_equal(fusion@gene_downstream@ensembl_id, correct_ensembl_id_downstream)

})

test_that("get_ensembl_ids gives error when it can't find ensembl id", {
  # Store the correct gene names temporarily
  name_upstream <- fusion@gene_upstream@name
  name_downstream <- fusion@gene_downstream@name

  # Set the geneNames[1] to a non-existant gene name
  fusion@gene_upstream@name <- "non-existant"
  fusion@gene_downstream@name <- "non-existant"

  # Check that we get a warning
  expect_error(get_ensembl_ids(fusion))

  # Restore original gene names
  fusion@gene_upstream@name <- name_upstream
  fusion@gene_downstream@name <- name_downstream
})

test_that("fusion_to_data_frame works as expected", {
  data_frame_fusion <- fusion_to_data_frame(fusion)
  expect_is(data_frame_fusion, "data.frame")
  expect_equal(
    fusion@gene_upstream@name,
    as.character(data_frame_fusion$gene_upstream)
  )

  # Expect error when given something other than a Fusion object
  expect_error(fusion_to_data_frame(42))
})

test_that("down_shift works as expected", {
  ir <- IRanges::IRanges(
    start = c(13, 40, 100),
    end = c(20, 53, 110))
  gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = ir)
  # Downshift it and see the introns are removed:
  gr_downshifted <- down_shift(gr)
  expect_is(gr_downshifted, "GRanges")
  expect_equal(IRanges::width(gr), IRanges::width(gr_downshifted))
})
