context("Test importing fusioncatcher data into a list of Fusion objects")

fusioncatcher833ke <- system.file(
  "extdata",
  "fusioncatcher_833ke_final-list-candidate-fusion-genes.txt",
  package = "chimeraviz")
fusions <- import_fusioncatcher(fusioncatcher833ke, "hg38", 3)
fusion <- fusions[[1]]

test_that("import_fusioncatcher imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "CLEC6A")
  expect_equal(fusion@gene_upstream@chromosome, "chr12")
})

test_that("import_fusioncatcher fails to import data", {
  expect_error(import_fusioncatcher("Non-existent-file", "hg38"))
  expect_error(import_fusioncatcher(fusioncatcher833ke, "Invalid-genome"))
  expect_error(
    import_fusioncatcher(fusioncatcher833ke, "hg38", "invalid limit"))
})

test_that("import_fusioncatcher limits work", {
  expect_equal(length(import_fusioncatcher(fusioncatcher833ke, "hg38", 2)), 2)
})
