context("Test importing SQUID data into a list of Fusion objects")

squidfile <- system.file(
  "extdata",
  "squid_hcc1954_sv.txt",
  package = "chimeraviz")
fusions <- import_squid(squidfile, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_squid imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@chromosome, "17")
})

test_that("import_squid fails to import data", {
  expect_error(import_squid("Non-existent-file", "hg18"))
  expect_error(import_squid(squidfile, "Invalid-genome"))
  expect_error(import_squid(squidfile, "hg19", "invalid limit"))
})

test_that("import_squid limits work", {
  expect_equal(length(import_squid(squidfile, "hg19", 2)), 2)
})
