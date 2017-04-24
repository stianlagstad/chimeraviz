#' An S4 class to represent a gene partner in a fusion
#'
#' The PartnerGene class represents one of the genes in a fusion event.
#'
#' @slot name Character containing name of the gene.
#' @slot ensemblId Character containing ensembl id for the gene.
#' @slot chromosome Character containing chromosome name.
#' @slot breakpoint Numeric containing the fusion breakpoint.
#' @slot strand Character containing gene strand.
#' @slot junctionSequence Biostrings::DNAString containing the sequence right before/after the fusion breakpoint.
#' @slot transcripts GenomicRanges::GRangesList containing three GenomicRanges::Granges() objects, one for each "transcript type". The transcript types are: 1) Transcripts where the fusion breakpoint hits an exon boundary, 2) transcripts where the fusion breakpoint is within an exon, 3) transcripts where the fusion breakpoint is within an intron.
#'
#' @export
PartnerGene <- setClass("PartnerGene",
                        slots = list(
                          name = "character",
                          ensemblId = "character",
                          chromosome = "character",
                          breakpoint = "numeric",
                          strand = "character",
                          junctionSequence = "DNAString",
                          transcripts = "GRangesList"
                        ))

#' Show method for the PartnerGene class.
#'
#' @param object A PartnerGene object
#'
#' @return Shows information about a PartnerGene object.
setMethod("show",
          "PartnerGene",
          function(object) {
            print("PartnerGene object")
            print(paste("Name:", object@name))
            print(paste("ensemblId:", object@ensemblId))
            print(paste("Chromosome:", object@chromosome))
            print(paste("Strand:", object@strand))
            print(paste("Breakpoint:", object@breakpoint))
          }
)

#' An S4 class to represent a fusion event.
#'
#' The Fusion class represents a fusion event, holding data imported from a
#' fusion tool.
#'
#' @slot id A unique id representing a fusion event. For deFuse data this will
#' be the cluster id.
#' @slot fusionTool Name of the fusion tool.
#' @slot genomeVersion Name of the genome used to map reads.
#' @slot spanningReadsCount The number of spanning reads supporting the fusion.
#' @slot splitReadsCount The number of split reads supporting the fusion.
#' @slot fusionReadsAlignment A Gviz::AlignmentsTrack object holding the fusion
#' reads aligned to the fusion sequence.
#' @slot geneA A PartnerGene object holding information of the upstream fusion
#' partner gene.
#' @slot geneB A PartnerGene object holding information of the downstream fusion
#' partner gene.
#' @slot inframe A logical value indicating whether or not the downstream fusion
#' partner gene is inframe or not. Not all fusion-finders report this.
#' @slot fusionToolSpecificData A list that will hold fields of importance for a
#' specific fusion finder. This field is used because many fusion-finders report
#' important values that are hard to fit into a standardized format. Examples of
#' values that are added to this list is probability from deFuse and EricScore
#' from EricScript.
#'
#' @export
Fusion <- setClass("Fusion",
                   slots = list(
                     id = "character",
                     fusionTool = "character",
                     genomeVersion = "character",
                     spanningReadsCount = "numeric",
                     splitReadsCount = "numeric",
                     fusionReadsAlignment = "AlignmentsTrack",
                     geneA = "PartnerGene",
                     geneB = "PartnerGene",
                     inframe = "logical",
                     fusionToolSpecificData = "list"
                   ))

#' Show method for the Fusion class.
#'
#' @param object A Fusion object
#'
#' @return Shows information about a Fusion object.
setMethod("show",
          "Fusion",
          function(object) {
            print("Fusion object")
            print(paste("id:", object@id))
            print(paste("Fusion tool:", object@fusionTool))
            print(paste("Genome version:", object@genomeVersion))
            print(paste("Gene names: ",
                        paste(object@geneA@name,
                              object@geneB@name,
                              sep = "-"),
                        sep = ""))
            print(paste("Chromosomes: ",
                        paste(object@geneA@chromosome,
                              object@geneB@chromosome,
                              sep = "-"),
                        sep = ""))
            print(paste("Strands: ",
                        paste(object@geneA@strand,
                              object@geneB@strand,
                              sep = ","),
                        sep = ""))
            print(paste("In-frame?: ",
                        object@inframe,
                        sep = ""))
          }
)

#' Get the upstream fusion partner gene
#'
#' This getter retrieves the upstream PartnerGene object.
#'
#' @param x The Fusion object you wish to retrieve the upstream PartnerGene object for.
#'
#' @return The upstream PartnerGene object.
#'
#' @rdname upstreamPartnerGene
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the upstream fusion partner gene
#' upstreamPartnerGene(fusion)
#'
#' @export
setGeneric(
  "upstreamPartnerGene",
  function(x)
    standardGeneric("upstreamPartnerGene")
)
#' @rdname upstreamPartnerGene
setMethod(
  "upstreamPartnerGene",
  "Fusion",
  function(x)
    return(x@geneA)
)

#' Get the downstream fusion partner gene
#'
#' This getter retrieves the downstream PartnerGene object.
#'
#' @param x The Fusion object you wish to retrieve the downstream PartnerGene object for.
#'
#' @return The downstream PartnerGene object.
#'
#' @rdname downstreamPartnerGene
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the downstream fusion partner gene
#' downstreamPartnerGene(fusion)
#'
#' @export
setGeneric(
  "downstreamPartnerGene",
  function(x)
    standardGeneric("downstreamPartnerGene")
)
#' @rdname downstreamPartnerGene
setMethod(
  "downstreamPartnerGene",
  "Fusion",
  function(x)
    return(x@geneB)
)

#' Set the upstream PartnerGene object of a Fusion object
#'
#' This sets the upstream PartnerGene object of a Fusion object
#'
#' @param object The Fusion object you wish to set a new upstream PartnerGene object for.
#' @param value The new PartnerGene object.
#'
#' @rdname upstreamPartnerGene
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Set the upstream PartnerGene object to be the same as the downstream PartnerGene object
#' upstreamPartnerGene(fusion) <- downstreamPartnerGene(fusion)
#'
#' @export
setGeneric(
  "upstreamPartnerGene<-",
  function(object,value) {
    standardGeneric("upstreamPartnerGene<-")
  }
)
#' @rdname upstreamPartnerGene
setReplaceMethod(
  f = "upstreamPartnerGene",
  signature = "Fusion",
  definition = function(object, value) {
    object@geneA <-value
    return (object)
  }
)

#' Set the downstream PartnerGene object of a Fusion object
#'
#' This sets the downstream PartnerGene object of a Fusion object
#'
#' @param object The Fusion object you wish to set a new downstream PartnerGene object for.
#' @param value The new PartnerGene object.
#'
#' @rdname downstreamPartnerGene
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Set the downstream PartnerGene object to be the same as the upstream PartnerGene object
#' downstreamPartnerGene(fusion) <- upstreamPartnerGene(fusion)
#'
#' @export
setGeneric(
  "downstreamPartnerGene<-",
  function(object,value) {
    standardGeneric("downstreamPartnerGene<-")
  }
)
#' @rdname downstreamPartnerGene
setReplaceMethod(
  f = "downstreamPartnerGene",
  signature = "Fusion",
  definition = function(object, value) {
    object@geneB <-value
    return (object)
  }
)

#' Get the Ensembl ID from a PartnerGene object
#'
#' This getter retrieves the Ensembl ID from a PartnerGene object
#'
#' @param x The PartnerGene object you wish to retrieve the Ensembl ID for.
#'
#' @return The upstream fusion partner gene Ensembl ID.
#'
#' @rdname partnerGeneEnsemblId
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the Ensembl ID from the upstream fusion partner gene
#' partnerGeneEnsemblId(upstreamPartnerGene(fusion))
#'
#' @export
setGeneric(
  "partnerGeneEnsemblId",
  function(x)
    standardGeneric("partnerGeneEnsemblId")
)
#' @rdname partnerGeneEnsemblId
setMethod(
  "partnerGeneEnsemblId",
  "PartnerGene",
  function(x)
    return(x@ensemblId)
)

#' Set the Ensembl ID of a PartnerGene object
#'
#' This sets the Ensembl ID of a PartnerGene object.
#'
#' @param object The PartnerGene object you wish to set a new Ensembl ID for.
#' @param value The new Ensembl ID.
#'
#' @rdname partnerGeneEnsemblId
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Set the downstream PartnerGene object to be the same as the upstream PartnerGene object
#' partnerGeneEnsemblId(upstreamPartnerGene(fusion)) <- "test"
#'
#' @export
setGeneric(
  "partnerGeneEnsemblId<-",
  function(object,value) {
    standardGeneric("partnerGeneEnsemblId<-")
  }
)
#' @rdname partnerGeneEnsemblId
setReplaceMethod(
  f = "partnerGeneEnsemblId",
  signature = "PartnerGene",
  definition = function(object, value) {
    object@ensemblId <-value
    return (object)
  }
)

#' Get the junction sequence from a PartnerGene object
#'
#' This getter retrieves the junction sequence from a PartnerGene object
#'
#' @param x The PartnerGene object you wish to retrieve the junction sequence for.
#'
#' @return The upstream fusion partner gene junction sequence.
#'
#' @rdname partnerGeneJunctionSequence
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the junction sequence from the upstream fusion partner gene
#' partnerGeneJunctionSequence(upstreamPartnerGene(fusion))
#'
#' @export
setGeneric(
  "partnerGeneJunctionSequence",
  function(x)
    standardGeneric("partnerGeneJunctionSequence")
)
#' @rdname partnerGeneJunctionSequence
setMethod(
  "partnerGeneJunctionSequence",
  "PartnerGene",
  function(x)
    return(x@junctionSequence)
)

#' Get the split reads count from a Fusion object
#'
#' This getter retrieves the split reads count from a Fusion object
#'
#' @param x The Fusion object you wish to retrieve the split reads count for.
#'
#' @return The Fusion split reads count.
#'
#' @rdname fusionSplitReadsCount
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the split reads count
#' fusionSplitReadsCount(fusion)
#'
#' @export
setGeneric(
  "fusionSplitReadsCount",
  function(x)
    standardGeneric("fusionSplitReadsCount")
)
#' @rdname fusionSplitReadsCount
setMethod(
  "fusionSplitReadsCount",
  "Fusion",
  function(x)
    return(x@splitReadsCount)
)

#' Get the spanning reads count from a Fusion object
#'
#' This getter retrieves the spanning reads count from a Fusion object
#'
#' @param x The Fusion object you wish to retrieve the spanning reads count for.
#'
#' @return The Fusion spanning reads count.
#'
#' @rdname fusionSpanningReadsCount
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the spanning reads count
#' fusionSpanningReadsCount(fusion)
#'
#' @export
setGeneric(
  "fusionSpanningReadsCount",
  function(x)
    standardGeneric("fusionSpanningReadsCount")
)
#' @rdname fusionSpanningReadsCount
setMethod(
  "fusionSpanningReadsCount",
  "Fusion",
  function(x)
    return(x@spanningReadsCount)
)
