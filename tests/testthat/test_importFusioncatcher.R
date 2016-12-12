context("Test importing fusioncatcher data into a list of Fusion objects")

fusioncatcher833ke <- system.file(
  "extdata",
  "fusioncatcher_833ke_final-list-candidate-fusion-genes.txt",
  package = "chimeraviz")
fusions <- importFusioncatcher(fusioncatcher833ke, "hg38", 3)
fusion <- fusions[[1]]

test_that("importFusioncatcher imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@geneA@name, "CLEC6A")
  expect_equal(fusion@geneA@chromosome, "chr12")
})

test_that("importFusioncatcher fails to import data", {
  expect_error(importFusioncatcher("Non-existent-file", "hg38"))
  expect_error(importFusioncatcher(fusioncatcher833ke, "Invalid-genome"))
  expect_error(importFusioncatcher(fusioncatcher833ke, "hg38", "invalid limit"))
})

test_that("importFusioncatcher limits work", {
  expect_equal(length(importFusioncatcher(fusioncatcher833ke, "hg38", 2)), 2)
})
