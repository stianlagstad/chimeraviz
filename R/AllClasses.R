#' An S4 class to represent a gene partner in a fusion
#'
#' The PartnerGene class represents one of the genes in a fusion event.
#'
#' @slot name Character containing name of the gene.
#' @slot ensembl_id Character containing ensembl id for the gene.
#' @slot chromosome Character containing chromosome name.
#' @slot breakpoint Numeric containing the fusion breakpoint.
#' @slot strand Character containing gene strand.
#' @slot junction_sequence Biostrings::DNAString containing the sequence right before/after the fusion breakpoint.
#' @slot transcripts GenomicRanges::GRangesList containing three GenomicRanges::Granges() objects, one for each "transcript type". The transcript types are: 1) Transcripts where the fusion breakpoint hits an exon boundary, 2) transcripts where the fusion breakpoint is within an exon, 3) transcripts where the fusion breakpoint is within an intron.
#'
#' @export
PartnerGene <- setClass(
  "PartnerGene",
  slots = c(
    name = "character",
    ensembl_id = "character",
    chromosome = "character",
    breakpoint = "numeric",
    strand = "character",
    junction_sequence = "DNAString",
    transcripts = "GRangesList"
  )
)

#' Show method for the PartnerGene class.
#'
#' @param object A PartnerGene object
#'
#' @return Shows information about a PartnerGene object.
setMethod(
  "show",
  "PartnerGene",
  function(object) {
    print("PartnerGene object")
    print(paste("Name:", object@name))
    print(paste("ensemblId:", object@ensembl_id))
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
#' @slot fusion_tool Name of the fusion tool.
#' @slot genome_version Name of the genome used to map reads.
#' @slot spanning_reads_count The number of spanning reads supporting the
#' fusion.
#' @slot split_reads_count The number of split reads supporting the fusion.
#' @slot fusion_reads_alignment A Gviz::AlignmentsTrack object holding the
#' fusion reads aligned to the fusion sequence.
#' @slot gene_upstream A PartnerGene object holding information of the upstream
#' fusion partner gene.
#' @slot gene_downstream A PartnerGene object holding information of the
#' downstream fusion partner gene.
#' @slot inframe A logical value indicating whether or not the downstream
#' fusion partner gene is inframe or not. Not all fusion-finders report this.
#' @slot fusion_tool_specific_data A list that will hold fields of importance
#' for a specific fusion finder. This field is used because many fusion-finders
#' report important values that are hard to fit into a standardized format.
#' Examples of values that are added to this list is probability from deFuse
#' and EricScore from EricScript.
#'
#' @export
Fusion <- setClass(
  "Fusion",
  slots = list(
    id = "character",
    fusion_tool = "character",
    genome_version = "character",
    spanning_reads_count = "numeric",
    split_reads_count = "numeric",
    fusion_reads_alignment = "AlignmentsTrack",
    gene_upstream = "PartnerGene",
    gene_downstream = "PartnerGene",
    inframe = "logical",
    fusion_tool_specific_data = "list"
  )
)

#' Show method for the Fusion class.
#'
#' @param object A Fusion object
#'
#' @return Shows information about a Fusion object.
setMethod(
  "show",
  "Fusion",
  function(object) {
    print("Fusion object")
    print(paste("id:", object@id))
    print(paste("Fusion tool:", object@fusion_tool))
    print(paste("Genome version:", object@genome_version))
    print(
      paste(
        "Gene names: ",
        paste(
          object@gene_upstream@name,
          object@gene_downstream@name,
          sep = "-"
        ),
        sep = ""
      )
    )
    print(
      paste(
        "Chromosomes: ",
        paste(
          object@gene_upstream@chromosome,
          object@gene_downstream@chromosome,
          sep = "-"
        ),
        sep = ""
      )
    )
    print(
      paste(
        "Strands: ",
        paste(
          object@gene_upstream@strand,
          object@gene_downstream@strand,
          sep = ","
        ),
        sep = ""
      )
    )
    print(
      paste(
        "In-frame?: ",
        object@inframe,
        sep = ""
      )
    )
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
#' @rdname upstream_partner_gene
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the upstream fusion partner gene
#' upstream_partner_gene(fusion)
#'
#' @export
setGeneric(
  "upstream_partner_gene",
  function(x)
    standardGeneric("upstream_partner_gene")
)
#' @rdname upstream_partner_gene
setMethod(
  "upstream_partner_gene",
  "Fusion",
  function(x)
    return(x@gene_upstream)
)

#' Get the downstream fusion partner gene
#'
#' This getter retrieves the downstream PartnerGene object.
#'
#' @param x The Fusion object you wish to retrieve the downstream PartnerGene object for.
#'
#' @return The downstream PartnerGene object.
#'
#' @rdname downstream_partner_gene
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the downstream fusion partner gene
#' downstream_partner_gene(fusion)
#'
#' @export
setGeneric(
  "downstream_partner_gene",
  function(x)
    standardGeneric("downstream_partner_gene")
)
#' @rdname downstream_partner_gene
setMethod(
  "downstream_partner_gene",
  "Fusion",
  function(x)
    return(x@gene_downstream)
)

#' Set the upstream PartnerGene object of a Fusion object
#'
#' This sets the upstream PartnerGene object of a Fusion object
#'
#' @param object The Fusion object you wish to set a new upstream PartnerGene object for.
#' @param value The new PartnerGene object.
#'
#' @rdname upstream_partner_gene
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Set the upstream PartnerGene object to be the same as the downstream PartnerGene object
#' upstream_partner_gene(fusion) <- downstream_partner_gene(fusion)
#'
#' @export
setGeneric(
  "upstream_partner_gene<-",
  function(object, value) {
    standardGeneric("upstream_partner_gene<-")
  }
)
#' @rdname upstream_partner_gene
setReplaceMethod(
  f = "upstream_partner_gene",
  signature = "Fusion",
  definition = function(object, value) {
    object@gene_upstream <- value
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
#' @rdname downstream_partner_gene
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Set the downstream PartnerGene object to be the same as the upstream PartnerGene object
#' downstream_partner_gene(fusion) <- upstream_partner_gene(fusion)
#'
#' @export
setGeneric(
  "downstream_partner_gene<-",
  function(object, value) {
    standardGeneric("downstream_partner_gene<-")
  }
)
#' @rdname downstream_partner_gene
setReplaceMethod(
  f = "downstream_partner_gene",
  signature = "Fusion",
  definition = function(object, value) {
    object@gene_downstream <- value
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
#' @rdname partner_gene_ensembl_id
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the Ensembl ID from the upstream fusion partner gene
#' partner_gene_ensembl_id(upstream_partner_gene(fusion))
#'
#' @export
setGeneric(
  "partner_gene_ensembl_id",
  function(x)
    standardGeneric("partner_gene_ensembl_id")
)
#' @rdname partner_gene_ensembl_id
setMethod(
  "partner_gene_ensembl_id",
  "PartnerGene",
  function(x)
    return(x@ensembl_id)
)

#' Set the Ensembl ID of a PartnerGene object
#'
#' This sets the Ensembl ID of a PartnerGene object.
#'
#' @param object The PartnerGene object you wish to set a new Ensembl ID for.
#' @param value The new Ensembl ID.
#'
#' @rdname partner_gene_ensembl_id
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Set the downstream PartnerGene object to be the same as the upstream PartnerGene object
#' partner_gene_ensembl_id(upstream_partner_gene(fusion)) <- "test"
#'
#' @export
setGeneric(
  "partner_gene_ensembl_id<-",
  function(object, value) {
    standardGeneric("partner_gene_ensembl_id<-")
  }
)
#' @rdname partner_gene_ensembl_id
setReplaceMethod(
  f = "partner_gene_ensembl_id",
  signature = "PartnerGene",
  definition = function(object, value) {
    object@ensembl_id <- value
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
#' @rdname partner_gene_junction_sequence
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the junction sequence from the upstream fusion partner gene
#' partner_gene_junction_sequence(upstream_partner_gene(fusion))
#'
#' @export
setGeneric(
  "partner_gene_junction_sequence",
  function(x)
    standardGeneric("partner_gene_junction_sequence")
)
#' @rdname partner_gene_junction_sequence
setMethod(
  "partner_gene_junction_sequence",
  "PartnerGene",
  function(x)
    return(x@junction_sequence)
)

#' Get the split reads count from a Fusion object
#'
#' This getter retrieves the split reads count from a Fusion object
#'
#' @param x The Fusion object you wish to retrieve the split reads count for.
#'
#' @return The Fusion split reads count.
#'
#' @rdname fusion_split_reads_count
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the split reads count
#' fusion_split_reads_count(fusion)
#'
#' @export
setGeneric(
  "fusion_split_reads_count",
  function(x)
    standardGeneric("fusion_split_reads_count")
)
#' @rdname fusion_split_reads_count
setMethod(
  "fusion_split_reads_count",
  "Fusion",
  function(x)
    return(x@split_reads_count)
)

#' Get the spanning reads count from a Fusion object
#'
#' This getter retrieves the spanning reads count from a Fusion object
#'
#' @param x The Fusion object you wish to retrieve the spanning reads count for.
#'
#' @return The Fusion spanning reads count.
#'
#' @rdname fusion_spanning_reads_count
#'
#' @examples
#' # Load data
#' defuseData <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuseData, "hg19", 1)
#' fusion <- fusions[[1]]
#' # Get the spanning reads count
#' fusion_spanning_reads_count(fusion)
#'
#' @export
setGeneric(
  "fusion_spanning_reads_count",
  function(x)
    standardGeneric("fusion_spanning_reads_count")
)
#' @rdname fusion_spanning_reads_count
setMethod(
  "fusion_spanning_reads_count",
  "Fusion",
  function(x)
    return(x@spanning_reads_count)
)
