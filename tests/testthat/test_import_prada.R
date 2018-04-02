context("Test importing PRADA data into a list of Fusion objects")

prada_data <- system.file(
  "extdata",
  "PRADA.acc.fusion.fq.TAF.tsv",
  package = "chimeraviz")
fusions <- import_prada(prada_data, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_prada imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "TBCD")
  expect_equal(fusion@gene_upstream@chromosome, "chr17")
  expect_equal(fusion@gene_upstream@strand, "+")
})

test_that("import_prada fails to import data", {
  expect_error(import_prada("Non-existent-file", "hg18"))
  expect_error(import_prada(prada_data, "Invalid-genome"))
  expect_error(import_prada(prada_data, "hg19", "invalid limit"))
})

test_that("import_prada limits work", {
  expect_equal(length(import_prada(prada_data, "hg19", 2)), 2)
})
