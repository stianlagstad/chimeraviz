context(
  paste0(
    "Test creating link data that RCircos will accept from a list of Fusion ",
    "object"
  )
)

defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 3)
link_data <- .fusions_to_link_data(fusions)

test_that(".fusions_to_link_data works", {
  expect_equal(length(link_data$chromosome), length(fusions))
  expect_equal(
    as.character(link_data$chromosome[1]),
    fusions[[1]]@gene_upstream@chromosome
  )
})
