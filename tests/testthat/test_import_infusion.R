context("Test importing InFusion data into a list of Fusion objects")

infusion_data <- system.file(
  "extdata",
  "infusion_fusions.txt",
  package = "chimeraviz")
fusions <- import_infusion(infusion_data, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_infusion imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "283813")
  expect_equal(fusion@gene_upstream@name, "TMPRSS2")
  expect_equal(fusion@gene_upstream@chromosome, "chr21")
})

test_that("import_infusion fails to import data", {
  expect_error(import_infusion("Non-existent-file", "hg18"))
  expect_error(import_infusion(infusion_data, "Invalid-genome"))
  expect_error(import_infusion(infusion_data, "hg19", "invalid limit"))
})

test_that("import_infusion limits work", {
  expect_equal(length(import_infusion(infusion_data, "hg19", 2)), 2)
})
