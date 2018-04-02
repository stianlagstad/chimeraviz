context(
  paste0(
    "Test creating gene label data that RCircos will accept from a list of ",
    "Fusion object"
  )
)

defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 3)
label_data <- .fusions_to_gene_label_data(fusions)

test_that(".fusions_to_gene_label_data works", {
  expect_equal(length(label_data$chromosome), length(fusions) * 2)
  expect_equal(
    as.character(label_data$chromosome[1]),
    fusions[[1]]@gene_upstream@chromosome
  )
})
