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
  grid::pushViewport(grid::viewport(layout = nf))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  Gviz::plotTracks(
    fusion@fusionReadsAlignment,
    from = fusion@geneA@breakpoint-1000,
    to = fusion@geneB@breakpoint+1000,
    chromosome = fusion@fusionReadsAlignment@chromosome,
    type = "pileup",
    add = FALSE)
  grid::popViewport(2)
}                                       #plotFusionReadsSimple
