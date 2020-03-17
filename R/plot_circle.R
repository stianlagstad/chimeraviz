#' Scale a vector of numeric values to an interval.
#'
#' This function takes a vector of numeric values as well as an interval
#' [new_min, new_max] that the numeric values will be scaled (normalized) to.
#'
#' @param the_list A vector of numeric values.
#' @param new_min Minimum value for the new interval.
#' @param new_max Maximum value for the new interval.
#'
#' @return A data frame with fusion link data compatible with RCircos::RCircos.Link.Plot()
#'
#' # @examples # Apparently examples shouldn't be set on private functions
#' list012 <- c(0,1,2)
#' .scale_list_to_interval(list012, 1, 3)
#' # [1] 1 2 3
#'
#' @name chimeraviz-internals-scaleListToInterval
.scale_list_to_interval <- function(the_list, new_min, new_max) {
if (length(the_list) <= 1) {
    stop(paste("Invalid list. Using this function with less than two values",
                "makes no sense."))
  }
  (new_max - new_min) * (the_list - min(the_list)) /
    (max(the_list) - min(the_list)) + new_min
}

#' Create link data for RCircos from the given fusions.
#'
#' This function takes a list of Fusion objects and creates a data frame in the
#' format that RCircos::RCircos.Link.Plot() expects for link data.
#'
#' @param fusion_list A list of Fusion objects.
#' @param min_link_width The minimum link line width. Default = 1
#' @param max_link_widt The maximum link line width. Default = 10
#'
#' @return A data frame with fusion link data compatible with RCircos::RCircos.Link.Plot()
#'
#' # @examples # Apparently examples shouldn't be set on private functions
#' defuse833ke <- system.file("extdata", "defuse_833ke_results.filtered.tsv", package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 3)
#' linkData <- chimeraviz::.fusions_to_link_data(fusions)
#' # This linkData can be used with RCircos::RCircos.Link.Plot()
#'
#' @name chimeraviz-internals-fusions_to_link_data
.fusions_to_link_data <- function(
  fusion_list,
  min_link_width = 1,
  max_link_widt = 10
) {
  chromosome <- vector(mode = "character", length = length(fusion_list))
  chrom_start <- vector(mode = "numeric", length = length(fusion_list))
  chrom_end <- vector(mode = "numeric", length = length(fusion_list))

  chromosome.1 <- vector(mode = "character", length = length(fusion_list))
  chrom_start.1 <- vector(mode = "numeric", length = length(fusion_list))
  chrom_end.1 <- vector(mode = "numeric", length = length(fusion_list))

  link_width <- vector(mode = "numeric", length = length(fusion_list))

  for (i in seq_along(fusion_list)) {
    fusion <- fusion_list[[i]]

    chromosome[[i]] <- fusion@gene_upstream@chromosome
    chrom_start[[i]] <- fusion@gene_upstream@breakpoint
    # This value shouldn't matter:
    chrom_end[[i]] <- fusion@gene_upstream@breakpoint + 1

    chromosome.1[[i]] <- fusion@gene_downstream@chromosome
    chrom_start.1[[i]] <- fusion@gene_downstream@breakpoint
    # This value shouldn't matter:
    chrom_end.1[[i]] <- fusion@gene_downstream@breakpoint + 1

    # Set link width = number of spanning+split reads supporting fusion
    link_width[[i]] <- fusion@spanning_reads_count + fusion@split_reads_count
  }

  # Normalize all link width values to the interval
  # [min_link_width, max_link_width]
  if (length(link_width) > 1) {
    link_width <-
      .scale_list_to_interval(link_width, min_link_width, max_link_widt)
  } else {
    link_width[[1]] <- max_link_widt
  }

  data.frame(chromosome,
             chrom_start,
             chrom_end,
             chromosome.1,
             chrom_start.1,
             chrom_end.1,
             link_width)
}

#' Create gene label data for RCircos from the given fusions.
#'
#' This function takes a list of Fusion objects and creates a data frame in the
#' format that RCircos.Gene.Name.Plot() expects for gene label data.
#'
#' @param fusion_list A list of Fusion objects.
#'
#' @return A data frame with fusion gene label data compatible with RCircos.Gene.Name.Plot()
#'
#' # @examples # Apparently examples shouldn't be set on private functions
#' defuse833ke <- system.file("extdata", "defuse_833ke_results.filtered.tsv", package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 3)
#' labelData <- chimeraviz::.fusions_to_gene_label_data(fusions)
#' # This labelData can be used with RCircos.Gene.Connector.Plot() and RCircos.Gene.Name.Plot()
#'
#' @name chimeraviz-internals-fusions_to_gene_label_data
.fusions_to_gene_label_data <- function(fusion_list) {
  original_length <- length(fusion_list)
  new_length <- original_length * 2

  chromosome <- vector(mode = "character", length = new_length)
  chrom_start <- vector(mode = "numeric", length = new_length)
  chrom_end <- vector(mode = "numeric", length = new_length)
  gene <- vector(mode = "character", length = new_length)

  for (i in seq_along(fusion_list)) {
    fusion <- fusion_list[[i]]

    # We are building the list of gene names for both partner genes at once, so
    # let's put the gene1 genes in the first half of the list..
    chromosome[[i]] <- fusion@gene_upstream@chromosome
    chrom_start[[i]] <- fusion@gene_upstream@breakpoint
    # This value shouldn't matter:
    chrom_end[[i]] <- fusion@gene_upstream@breakpoint + 1
    gene[[i]] <- fusion@gene_upstream@name

    # and put the gene2 genes in the other half.
    chromosome[[i + original_length]] <- fusion@gene_downstream@chromosome
    chrom_start[[i + original_length]] <- fusion@gene_downstream@breakpoint
    # This value shouldn't matter:
    chrom_end[[i + original_length]] <- fusion@gene_downstream@breakpoint + 1
    gene[[i + original_length]] <- fusion@gene_downstream@name
  }

  data.frame(chromosome,
             chrom_start,
             chrom_end,
             gene)
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
#' @param fusion_list A list of Fusion objects.
#'
#' @return Creates a circle plot.
#'
#' @examples
#' defuse833ke <- system.file(
#'   "extdata",
#'   "defuse_833ke_results.filtered.tsv",
#'   package="chimeraviz")
#' fusions <- import_defuse(defuse833ke, "hg19", 3)
#' # Temporary file to store the plot
#' pngFilename <- tempfile(
#'   pattern = "circlePlot",
#'   fileext = ".png",
#'   tmpdir = tempdir())
#' # Open device
#' png(pngFilename, width = 1000, height = 750)
#' # Plot!
#' plot_circle(fusions)
#' # Close device
#' dev.off()
#'
#' @export
plot_circle <- function(fusion_list) {

  .validate_plot_circle_params(fusion_list)

  # Read cytoband information depending on genome version
  if (fusion_list[[1]]@genome_version == "hg19") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.HG19.Human.CytoBandIdeogram.txt",
      package = "chimeraviz")
  } else if (fusion_list[[1]]@genome_version == "hg38") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.HG38.Human.CytoBandIdeogram.txt",
      package = "chimeraviz")
  } else if (fusion_list[[1]]@genome_version == "mm10") {
    cytoband_file <- system.file(
      "extdata",
      "UCSC.MM10.Mus.musculus.CytoBandIdeogram.txt",
      package = "chimeraviz"
    )
  } else {
    stop("Invalid genome version.")
  }
  cytoband <- utils::read.table(cytoband_file)
  # Set names to what RCircos expects
  names(cytoband) <- c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain")

  # We need the RCircos.Env object in the global namespace. Why? Because RCircos
  # is weird that way.
  assign("RCircos.Env", RCircos::RCircos.Env, .GlobalEnv)

  # Sort the cytoband data so that the chromosomes appear in order
  cytoband <- RCircos.Sort.Genomic.Data( # Exclude Linting
    genomic.data = cytoband,
    is.ideo = TRUE
  )

  # Initialize components
  cyto.info <- cytoband
  chr.exclude <- NULL
  tracks.inside <- 3
  tracks.outside <- 0
  RCircos::RCircos.Set.Core.Components(
    cyto.info,
    chr.exclude,
    tracks.inside,
    tracks.outside
  )

  # Open a new window for plotting
  RCircos::RCircos.Set.Plot.Area()

  # Draw chromosome ideogram
  RCircos::RCircos.Chromosome.Ideogram.Plot()

  # Create gene label data in the format RCircos requires
  gene_label_data <- .fusions_to_gene_label_data(fusion_list)

  # Draw gene names
  name_col <- 4;
  side <- "in";
  track_num <- 1;
  # Plot connectors
  RCircos.Gene.Connector.Plot( # Exclude Linting
    gene_label_data,
    track_num,
    side
  );
  track_num <- 2;
  # Plot gene names
  RCircos.Gene.Name.Plot( # Exclude Linting
    gene_label_data,
    name_col,
    track_num,
    side
  );

  # Create link data in the format RCircos requires
  link_data <- .fusions_to_link_data(fusion_list)
  # Make sure the ordering is correct.
  # Ref https://github.com/stianlagstad/chimeraviz/issues/52
  multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE) {
    do.call(
      order,
      c(
        lapply(
          list(...),
          function(l) {
            if (is.character(l)) {
              factor(
                l,
                levels = gtools::mixedsort(unique(l))
              )
            } else {
              l
            }
          }
        ),
        list(na.last = na.last, decreasing = decreasing)
      )
    )
  }
  ordered_link_width <- link_data[
    multi.mixedorder(
      as.character(link_data$chromosome),
      link_data$chrom_start
    ),
  ]$link_width

  # Draw link data
  # Which track?
  track_num <- 3
  # Plot links
  RCircos::RCircos.Link.Plot(
    link.data = link_data,
    track.num = track_num,
    by.chromosome = TRUE,
    start.pos = NULL,
    genomic.columns = 3,
    is.sorted = FALSE,
    lineWidth = ordered_link_width)

  # When done, remove the RCircos.Env element from the global namespace
  remove("RCircos.Env", envir = .GlobalEnv)
}

.validate_plot_circle_params <- function(
  fusion_list
) {
  # Establish a new 'ArgCheck' object
  argument_checker <- ArgumentCheck::newArgCheck()

  # Check input type and length
  if (
    !is.list(fusion_list) ||
      length(fusion_list) == 0
  ) {
    ArgumentCheck::addError(
      msg = "'fusion_list' must be a list of fusion objects",
      argcheck = argument_checker
    )
  }
  # Check that all the items in the input list are fusion objects
  for (i in seq_along(fusion_list)) {
    if (class(fusion_list[[i]]) != "Fusion") {
      ArgumentCheck::addError(
        msg = "'fusionList' must be a list of fusion objects",
        argcheck = argument_checker
      )
      break
    }
  }

  # Validate that the first fusion object has a valid genome
  genome_of_first_fusion <- fusion_list[[1]]@genome_version
  if (!genome_of_first_fusion %in% list("hg19", "hg38", "mm10")) {
    ArgumentCheck::addError(
      msg = paste0(
        "Invalid input. genomeVersion must be either \"hg19\", \"hg38\" or ",
        "\"mm10\"."
      ),
      argcheck = argument_checker
    )
  }
  # Validate that all the fusions have the same genome
  for (i in seq_along(fusion_list)) {
    if (!fusion_list[[i]]@genome_version == genome_of_first_fusion) {
      ArgumentCheck::addError(
        msg = paste0(
          "All Fusion objects in the fusion_list must have the same genome."
        ),
        argcheck = argument_checker
      )
    }
  }

  # Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(argument_checker)
}
