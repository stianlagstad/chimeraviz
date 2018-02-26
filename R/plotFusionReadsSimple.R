#' Simple versin of plotFusionReads
#'
#' 
#'
#' @export
plotFusionReadsSimple <- function(fusion, leftFlank=10000, rightFlank=leftFlank) {
  Gviz::displayPars(fusion@fusionReadsAlignment) <- list(
    showTitle = FALSE,
    showMismatches = FALSE # Show mismatched reads?
  )

  nf <- grid::grid.layout(
    nrow = 1,
    ncol = 1,
    heights = grid::unit(c(6), "null"),
    widths = grid::unit(c(1), "null"))
  
  grid::grid.newpage()
  Gviz::plotTracks(
    fusion@fusionReadsAlignment,
    from = fusion@geneA@breakpoint-leftFlank,
    to = fusion@geneB@breakpoint+rightFlank,
    chromosome = fusion@fusionReadsAlignment@chromosome,
    type = "pileup",
    add = TRUE)
}                                       #plotFusionReadsSimple
