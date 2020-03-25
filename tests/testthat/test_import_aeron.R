context("Test importing Aeron data into a list of Fusion objects")

aeronfusionsupportfile <- system.file(
 "extdata",
 "aeron_fusion_support.txt",
 package="chimeraviz")
aeronfusiontranscriptfile <- system.file(
 "extdata",
 "aeron_fusion_transcripts.fa",
 package="chimeraviz")
fusions <- import_aeron(
  aeronfusionsupportfile,
  aeronfusiontranscriptfile,
  "hg19",
  3)
fusion <- fusions[[1]]

test_that("import_aeron imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "115")
  expect_equal(fusion@gene_upstream@name, NA_character_)
  expect_equal(fusion@gene_upstream@chromosome, NA_character_)
})

test_that("import_aeron fails to import data", {
  expect_error(import_aeron("Non-existent-file", "hg18"))
  expect_error(
    import_aeron(
      aeronfusionsupportfile,
      aeronfusiontranscriptfile,
      "Invalid-genome"
    )
  )
  expect_error(
    import_aeron(
      aeronfusionsupportfile,
      aeronfusiontranscriptfile,
      "hg19",
      "invalid limit"
    )
  )
})

test_that("import_aeron limits work", {
  expect_equal(
    length(
      import_aeron(
        aeronfusionsupportfile,
        aeronfusiontranscriptfile,
        "hg19",
        2
      )
    ),
    2
  )
})

test_that(".aeron_parse_data works", {
  expect_equal(length(.aeron_parse_data("fusion_115_ENSG00000212907.2_291bp_ENSG00000198886.2_1340bp_50reads")), 4)
})

test_that(".import_aeron_transcript_data works", {
  expect_equal(length(.import_aeron_transcript_data(aeronfusiontranscriptfile)), 224)
})
