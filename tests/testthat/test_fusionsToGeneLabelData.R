context("Test creating gene label data that RCircos will accept from a list of Fusion object")

defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 3)
labelData <- .fusionsToGeneLabelData(fusions)

test_that("fusionsToGeneLabelData works", {
  expect_equal(length(labelData$Chromosome), length(fusions) * 2)
  expect_equal(as.character(labelData$Chromosome[1]), fusions[[1]]@geneA@chromosome)
})
