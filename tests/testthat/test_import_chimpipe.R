context("Test importing ChimPipe data into a list of Fusion objects")

chimpipefile <- system.file(
  "extdata",
  "chimericJunctions_MCF-7.txt",
  package = "chimeraviz")
fusions <- import_chimpipe(chimpipefile, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_chimpipe imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "KRT86")
  expect_equal(fusion@gene_upstream@chromosome, "chr12")
})

test_that("import_chimpipe fails to import data", {
  expect_error(import_chimpipe("Non-existent-file", "hg18"))
  expect_error(import_chimpipe(chimpipefile, "Invalid-genome"))
  expect_error(import_chimpipe(chimpipefile, "hg19", "invalid limit"))
})

test_that("import_chimpipe limits work", {
  expect_equal(length(import_chimpipe(chimpipefile, "hg19", 2)), 2)
})

test_that(".chimpipe_parse_junccoord works", {
  expect_equal(length(.chimpipe_parse_junccoord("chr12_9392928_+:chr11_75115716_+")), 6)
})
