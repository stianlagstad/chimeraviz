context("Test importing EricScript data into a list of Fusion objects")

ericscript_data <- system.file(
  "extdata",
  "ericscript_SRR1657556.results.total.tsv",
  package = "chimeraviz")
fusions <- import_ericscript(ericscript_data, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_ericscript imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "PPIA")
  expect_equal(fusion@gene_upstream@chromosome, "chr7")
})

test_that("import_ericscript fails to import data", {
  expect_error(import_ericscript("Non-existent-file", "hg18"))
  expect_error(import_ericscript(ericscript_data, "Invalid-genome"))
  expect_error(import_ericscript(ericscript_data, "hg19", "invalid limit"))
})

test_that("import_ericscript limits work", {
  expect_equal(length(import_ericscript(ericscript_data, "hg19", 2)), 2)
})
