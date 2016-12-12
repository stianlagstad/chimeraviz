#' Find a specific fusion object in a list by id
#'
#' Helper function to retrieve the Fusion object with the given id.
#'
#' @param fusionList A list of Fusion objects.
#' @param id The id (e.g. the cluster_id from a deFuse run) we're looking for.
#'
#' @return A Fusion object.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 1)
#' fusion <- getFusionById(fusions, 5267)
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
getFusionById <- function(fusionList, id) {
  result <- Filter(function(fusion) fusion@id==id, fusionList)
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
#' Note: getFusionByGeneName(fusionList, "MT") will match both MT-ND5 and
#' MT-ND4.
#'
#' @param fusionList A list of Fusion objects.
#' @param geneName The gene name we're looking for.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 1)
#' length(getFusionByGeneName(fusions, "RCC1"))
#' # [1] 1
#'
#' @export
getFusionByGeneName <- function(fusionList, geneName) {
  result <- Filter(function(fusion) {
    (grepl(geneName, fusion@geneA@name) == 1 ||
       grepl(geneName, fusion@geneB@name) == 1)
  }, fusionList)
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
#' @param fusionList A list of Fusion objects.
#' @param chr The chromosome name we're looking for fusions in.
#'
#' @return A list of Fusion objects.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 1)
#' length(getFusionByChromosome(fusions, "chr1"))
#' # [1] 1
#'
#' @export
getFusionByChromosome <- function(fusionList, chr) {
  result <- Filter(function(fusion) {
    (fusion@geneA@chromosome == chr ||
      fusion@geneB@chromosome == chr)
  }, fusionList)
  if (length(result) > 0) {
    return(result)
  } else {
    warning("Could not find fusions involving that chromosome.")
    return(result)
  }
}
