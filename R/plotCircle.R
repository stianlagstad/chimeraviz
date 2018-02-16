#' Scale a vector of numeric values to an interval.
#'
#' This function takes a vector of numeric values as well as an interval
#' [newMin, newMax] that the numeric values will be scaled (normalized) to.
#'
#' @param theList A vector of numeric values.
#' @param newMin Minimum value for the new interval.
#' @param newMax Maximum value for the new interval.
#'
#' @return A data frame with fusion link data compatible with RCircos::RCircos.Link.Plot()
#'
#' # @examples # Apparently examples shouldn't be set on private functions
#' list012 <- c(0,1,2)
#' .scaleListToInterval(list012, 1, 3)
#' # [1] 1 2 3
#'
#' @name chimeraviz-internals-scaleListToInterval
.scaleListToInterval <- function(theList, newMin, newMax) {
  if (length(theList) <= 1) {
    stop(paste("Invalid list. Using this function with less than two values",
                "makes to sense."))
  }
  (newMax-newMin)*(theList-min(theList))/(max(theList)-min(theList)) + newMin
}

#' Create link data for RCircos from the given fusions.
#'
#' This function takes a list of Fusion objects and creates a data frame in the
#' format that RCircos::RCircos.Link.Plot() expects for link data.
#'
#' @param fusionList A list of Fusion objects.
#' @param minLinkWidth The minimum link line width. Default = 1
#' @param maxLinkWidt The maximum link line width. Default = 10
#'
#' @return A data frame with fusion link data compatible with RCircos::RCircos.Link.Plot()
#'
#' # @examples # Apparently examples shouldn't be set on private functions
#' defuse833ke <- system.file("extdata", "defuse_833ke_results.filtered.tsv", package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 3)
#' linkData <- chimeraviz::.fusionsToLinkData(fusions)
#' # This linkData can be used with RCircos::RCircos.Link.Plot()
#'
#' @name chimeraviz-internals-fusionsToLinkData
.fusionsToLinkData <- function(fusionList, minLinkWidth = 1, maxLinkWidt = 10) {
  Chromosome <- vector(mode = "character", length = length(fusionList))
  chromStart <- vector(mode = "numeric", length = length(fusionList))
  chromEnd <- vector(mode = "numeric", length = length(fusionList))

  Chromosome.1 <- vector(mode = "character", length = length(fusionList))
  chromStart.1 <- vector(mode = "numeric", length = length(fusionList))
  chromEnd.1 <- vector(mode = "numeric", length = length(fusionList))

  linkWidth <- vector(mode = "numeric", length = length(fusionList))

  for (i in 1:length(fusionList)) {
    fusion <- fusionList[[i]]

    Chromosome[[i]] <- fusion@geneA@chromosome
    chromStart[[i]] <- fusion@geneA@breakpoint
    chromEnd[[i]] <- fusion@geneA@breakpoint + 1 # This value shouldn't matter

    Chromosome.1[[i]] <- fusion@geneB@chromosome
    chromStart.1[[i]] <- fusion@geneB@breakpoint
    chromEnd.1[[i]] <- fusion@geneB@breakpoint + 1 # This value shouldn't matter

    # Set link width = number of spanning+split reads supporting fusion
    linkWidth[[i]] <- fusion@spanningReadsCount + fusion@splitReadsCount
  }

  # Normalize all link width values to the interval [minLinkWidth, maxLinkWidt]
  if (length(linkWidth) > 1) {
    linkWidth <- .scaleListToInterval(linkWidth, minLinkWidth, maxLinkWidt)
  } else {
    linkWidth[[1]] <- maxLinkWidt
  }

  data.frame(Chromosome,
             chromStart,
             chromEnd,
             Chromosome.1,
             chromStart.1,
             chromEnd.1,
             linkWidth)
}

#' Create gene label data for RCircos from the given fusions.
#'
#' This function takes a list of Fusion objects and creates a data frame in the
#' format that RCircos.Gene.Name.Plot() expects for gene label data.
#'
#' @param fusionList A list of Fusion objects.
#'
#' @return A data frame with fusion gene label data compatible with RCircos.Gene.Name.Plot()
#'
#' # @examples # Apparently examples shouldn't be set on private functions
#' defuse833ke <- system.file("extdata", "defuse_833ke_results.filtered.tsv", package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 3)
#' labelData <- chimeraviz::.fusionsToGeneLabelData(fusions)
#' # This labelData can be used with RCircos.Gene.Connector.Plot() and RCircos.Gene.Name.Plot()
#'
#' @name chimeraviz-internals-fusionsToGeneLabelData
.fusionsToGeneLabelData <- function(fusionList) {
  originalLength <- length(fusionList)
  newLength <- originalLength * 2

  Chromosome <- vector(mode = "character", length = newLength)
  chromStart <- vector(mode = "numeric", length = newLength)
  chromEnd <- vector(mode = "numeric", length = newLength)
  Gene <- vector(mode = "character", length = newLength)

  for (i in 1:length(fusionList)) {
    fusion <- fusionList[[i]]

    # We are building the list of gene names for both partner genes at once, so
    # let's put the gene1 genes in the first half of the list..
    Chromosome[[i]] <- fusion@geneA@chromosome
    chromStart[[i]] <- fusion@geneA@breakpoint
    chromEnd[[i]] <- fusion@geneA@breakpoint + 1 # This value shouldn't matter
    Gene[[i]] <- fusion@geneA@name

    # and put the gene2 genes in the other half.
    Chromosome[[i+originalLength]] <- fusion@geneB@chromosome
    chromStart[[i+originalLength]] <- fusion@geneB@breakpoint
    chromEnd[[i+originalLength]] <- fusion@geneB@breakpoint + 1 # This value shouldn't matter
    Gene[[i+originalLength]] <- fusion@geneB@name
  }

  data.frame(Chromosome,
             chromStart,
             chromEnd,
             Gene)
}

#' Create a circle plot of the given fusions.
#'
#' This function takes a list of Fusion objects and creates a circle plot
#' indicating which chromosomes the fusion genes in the list consists of.
#'
#' Note that only a limited number of gene names can be shown in the circle plot
#' due to the limited resolution of the plot. RCircos will automatically limit
#' the number of gene names shown if there are too many. Also note that fusions
#' involving mitochondrial DNA will not be shown in this plot.
#'
#' @param fusionList A list of Fusion objects.
#'
#' @return Creates a circle plot.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- importDefuse(defuse833ke, "hg19", 3)
#' # Temporary file to store the plot
#' pngFilename <- tempfile(
#'   pattern = "circlePlot",
#'   fileext = ".png",
#'   tmpdir = tempdir())
#' # Open device
#' png(pngFilename, width = 1000, height = 750)
#' # Plot!
#' plotCircle(fusions)
#' # Close device
#' dev.off()
#'
#' @export
plotCircle <- function(fusionList) {

  # Check that input is a list
  if (!is.list(fusionList)) {
    stop("Input must be a list of fusion objects")
  }

  # Check that first item in input list is a fusion object
  if (class(fusionList[[1]]) != "Fusion") {
    stop("items in fusionList must be an object of type Fusion")
  }

  if (!fusionList[[1]]@genomeVersion %in% list("hg19", "hg38", "mm10")) {
    stop("Invalid input. genomeVersion must be either \"hg19\", \"hg38\" or \"mm10\".")
  }

##   if (any(rapply(fusionList, function(fusion) fusion@geneA@chromosome == "chrM" || fusion@geneB@chromosome == "chrM"))) {
##     indexesMitochondrialGenes <-
##       rapply(fusionList, function(fusion) fusion@geneA@chromosome == "chrM" || fusion@geneB@chromosome == "chrM")
##     fusionList <- fusionList[!indexesMitochondrialGenes]
##     message(
##       paste0(
##         "Removing ",
##         length(fusionList[indexesMitochondrialGenes]),
##         " fusions involving mitochondrial genes as they cannot be plotted in the circle plot."
##         )
##       )
##   }

  # Read cytoband information depending on genome version
  if (fusionList[[1]]@genomeVersion == "hg19") {
    cytobandFile <- system.file(
      "extdata",
      "UCSC.HG19.Human.CytoBandIdeogram.txt",
      package="chimeraviz")
  } else if (fusionList[[1]]@genomeVersion == "hg38") {
    cytobandFile <- system.file(
      "extdata",
      "UCSC.HG38.Human.CytoBandIdeogram.txt",
      package="chimeraviz")
  } else if (fusionList[[1]]@genomeVersion == "mm10") {
    cytobandFile <- system.file(
      "extdata",
      "UCSC.MM10.Mus.musculus.CytoBandIdeogram.txt",
      package="chimeraviz"
    )
  }
  cytoband <- utils::read.table(cytobandFile)
  # Set names to what RCircos expects
  names(cytoband) <- c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain")

  # We need the RCircos.Env object in the global namespace. Why? Because RCircos
  # is weird that way.
  assign("RCircos.Env", RCircos::RCircos.Env, .GlobalEnv)

  # Sort the cytoband data so that the chromosomes appear in order
  cytoband <- RCircos.Sort.Genomic.Data(genomic.data = cytoband, is.ideo = TRUE)

  # Initialize components
  cyto.info <- cytoband
  chr.exclude <- NULL
  tracks.inside <- 3
  tracks.outside <- 0
  RCircos::RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

  # Open a new window for plotting
  RCircos::RCircos.Set.Plot.Area()

  # Draw chromosome ideogram
  RCircos::RCircos.Chromosome.Ideogram.Plot()

  # Create gene label data in the format RCircos requires
  geneLabelData <- .fusionsToGeneLabelData(fusionList)

  # Draw gene names
  name.col <- 4;
  side <- "in";
  track.num <- 1;
  # Plot connectors
  RCircos.Gene.Connector.Plot(geneLabelData, track.num, side);
  track.num <- 2;
  # Plot gene names
  RCircos.Gene.Name.Plot(geneLabelData, name.col,track.num, side);

  # Create link data in the format RCircos requires
  linkData <- .fusionsToLinkData(fusionList)

  # Draw link data
  # Which track?
  track.num <- 3
  # Plot links
  RCircos::RCircos.Link.Plot(
    link.data = linkData,
    track.num = track.num,
    by.chromosome = TRUE,
    start.pos = NULL,
    genomic.columns = 3,
    is.sorted = FALSE,
    lineWidth = linkData$linkWidth)

  # When done, remove the RCircos.Env element from the global namespace
  remove("RCircos.Env", envir = .GlobalEnv)
}

.validatePlotCircleParams <- function(
  fusionList
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()

  # Add an error if the input isn't a list of fusion objects
  if (!is.list(fusionList) || length(fusionList) == 0 || class(fusionList[[1]]) != "Fusion") {
    ArgumentCheck::addError(
      msg = "'fusionList' must be a list of fusion objects",
      argcheck = argument_checker
    )
  }

  if (!fusionList[[1]]@genomeVersion %in% list("hg19", "hg38")) {
    ArgumentCheck::addError(
      msg = paste0("The genomeVersion of the Fusion objects must be either",
                   "'hg19' or 'hg38'."),
      argcheck = argument_checker
    )
  }

  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)
}
