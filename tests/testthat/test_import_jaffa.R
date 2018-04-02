context("Test importing JAFFA data into a list of Fusion objects")

jaffa_data <- system.file(
  "extdata",
  "jaffa_results.csv",
  package = "chimeraviz")
fusions <- import_jaffa(jaffa_data, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_jaffa imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "BCAS4")
  expect_equal(fusion@gene_upstream@chromosome, "chr20")
  expect_equal(fusion@gene_upstream@strand, "*")
})

test_that("import_jaffa fails to import data", {
  expect_error(import_jaffa("Non-existent-file", "hg18"))
  expect_error(import_jaffa(jaffa_data, "Invalid-genome"))
  expect_error(import_jaffa(jaffa_data, "hg19", "invalid limit"))
})

test_that("import_jaffa limits work", {
  expect_equal(length(import_jaffa(jaffa_data, "hg19", 2)), 2)
})
