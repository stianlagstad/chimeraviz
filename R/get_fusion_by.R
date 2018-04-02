#' Find a specific fusion object in a list by id
#'
#' Helper function to retrieve the Fusion object with the given id.
#'
#' @param fusion_list A list of Fusion objects.
#' @param id The id (e.g. the cluster_id from a deFuse run) we're looking for.
#'
#' @return A Fusion object.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' fusion <- get_fusion_by_id(fusions, 5267)
#' # This should be the Fusion object:
#' fusion
#' # [1] "Fusion object"
#' # [1] "id: 5267"
#' # [1] "Fusion tool: defuse"
#' # [1] "Genome version: hg19"
#' # [1] "Gene names: RCC1-HENMT1"
#' # [1] "Chromosomes: chr1-chr1"
#' # [1] "Strands: +,-"
#'
#' @export
get_fusion_by_id <- function(fusion_list, id) {
  result <- Filter(function(fusion) fusion@id == id, fusion_list)
  if (length(result) == 1) {
    return(result[[1]])
  } else {
    stop("Could not find fusion.")
  }
}

#' Find fusions that includes the given gene.
#'
#' Helper function to retrieve the Fusion objects that has geneName as one of
#' the partner genes.
#'
#' Note: get_fusion_by_gene_name(fusionList, "MT") will match both MT-ND5 and
#' MT-ND4.
#'
#' @param fusion_list A list of Fusion objects.
#' @param gene_name The gene name we're looking for.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' length(get_fusion_by_gene_name(fusions, "RCC1"))
#' # [1] 1
#'
#' @export
get_fusion_by_gene_name <- function(fusion_list, gene_name) {
  result <- Filter(
    function(fusion) {
      (grepl(gene_name, fusion@gene_upstream@name) == 1 ||
         grepl(gene_name, fusion@gene_downstream@name) == 1)
      },
    fusion_list
  )
  if (length(result) > 0) {
    return(result)
  } else {
    warning("Could not find fusions involving that gene name.")
    return(result)
  }
}

#' Find fusions that involves genes in the given chromosome.
#'
#' Helper function to retrieve the Fusion objects that involves genes in the
#' given chromosome name.
#'
#' @param fusion_list A list of Fusion objects.
#' @param chr The chromosome name we're looking for fusions in.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 1)
#' length(get_fusion_by_chromosome(fusions, "chr1"))
#' # [1] 1
#'
#' @export
get_fusion_by_chromosome <- function(fusion_list, chr) {
  result <- Filter(
    function(fusion) {
      (fusion@gene_upstream@chromosome == chr ||
         fusion@gene_downstream@chromosome == chr)
      },
    fusion_list
  )
  if (length(result) > 0) {
    return(result)
  } else {
    warning("Could not find fusions involving that chromosome.")
    return(result)
  }
}
