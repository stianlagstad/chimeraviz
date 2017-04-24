context("Test importing defuse data into a list of Fusion objects")

defuse833ke <- system.file("extdata",
                           "defuse_833ke_results.filtered.tsv",
                           package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 3)
fusion <- fusions[[1]]

test_that("importDefuse imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "5267")
  expect_equal(fusion@geneA@name, "RCC1")
  expect_equal(fusion@geneA@chromosome, "chr1")
})

test_that("importDefuse fails to import data", {
  expect_error(importDefuse("Non-existent-file", "hg18"))
  expect_error(importDefuse(defuse833ke, "Invalid-genome"))
  expect_error(importDefuse(defuse833ke, "hg19", "invalid limit"))
})

test_that("importDefuse limits work", {
  expect_equal(length(importDefuse(defuse833ke, "hg19", 2)), 2)
})
