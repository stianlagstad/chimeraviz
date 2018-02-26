#' Simple versin of plotFusionReads
#'
#' 
#'
#' @export
plotFusionReadsSimple <- function(fusion) {
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
    from = fusion@geneA@breakpoint-10000,
    to = fusion@geneB@breakpoint+10000,
    chromosome = fusion@fusionReadsAlignment@chromosome,
    type = "pileup",
    add = TRUE)
}                                       #plotFusionReadsSimple
