context("Test creating link data that RCircos will accept from a list of Fusion object")

defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 3)
linkData <- .fusionsToLinkData(fusions)

test_that("fusionsToLinkData works", {
  expect_equal(length(linkData$Chromosome), length(fusions))
  expect_equal(as.character(linkData$Chromosome[1]), fusions[[1]]@geneA@chromosome)
})
