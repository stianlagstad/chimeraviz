context("Test importing defuse data into a list of Fusion objects")

defuse833ke <- system.file("extdata",
                           "defuse_833ke_results.filtered.tsv",
                           package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_defuse imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "5267")
  expect_equal(fusion@gene_upstream@name, "RCC1")
  expect_equal(fusion@gene_upstream@chromosome, "chr1")
})

test_that("import_defuse fails to import data", {
  expect_error(import_defuse("Non-existent-file", "hg18"))
  expect_error(import_defuse(defuse833ke, "Invalid-genome"))
  expect_error(import_defuse(defuse833ke, "hg19", "invalid limit"))
})

test_that("import_defuse limits work", {
  expect_equal(length(import_defuse(defuse833ke, "hg19", 2)), 2)
})
