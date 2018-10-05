context("Test importing oncofuse data into a list of Fusion objects")

oncofuse_file <- system.file(
  "extdata",
  "oncofuse.outfile",
  package = "chimeraviz")
fusions <- import_oncofuse(oncofuse_file, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_oncofuse imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "PGC")
  expect_equal(fusion@gene_upstream@chromosome, "chr6")
})

test_that("import_oncofuse fails to import data", {
  expect_error(import_oncofuse("Non-existent-file", "hg18"))
  expect_error(import_oncofuse(oncofuse_file, "Invalid-genome"))
  expect_error(import_oncofuse(oncofuse_file, "hg19", "invalid limit"))
})

test_that("import_oncofuse limits work", {
  expect_equal(length(import_oncofuse(oncofuse_file, "hg19", 2)), 2)
})
