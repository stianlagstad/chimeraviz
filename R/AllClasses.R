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
