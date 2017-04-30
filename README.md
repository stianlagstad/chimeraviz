# chimeraviz

`chimeraviz` is an R package that automates the creation of chimeric RNA visualizations. See the [package vignette](https://bioconductor.org/packages/release/bioc/vignettes/chimeraviz/inst/doc/chimeraviz-vignette.html) for more information.

# Installing chimeraviz

`chimeraviz` is a Bioconductor package, and is most easily installed via Bioconductor. See instructions on how to do that [here](https://bioconductor.org/packages/release/bioc/html/chimeraviz.html).

If you would like to build the package yourself, first install these CRAN dependencies:

```
# CRAN
cranPackages <- c(
  "devtools",     # For installing chimeraviz
  "graphics",     # For plotting
  "DT",           # For creating a fusion report
  "rmarkdown",    # For creating a fusion report
  "readr",        # For reading fusion-finder result files
  "grid",         # For plotting
  "RColorBrewer", # For picking colors in plotting
  "RCircos",      # For creating overview plot
  "plyr",         # For converting fusion data to data.frame format
  "dplyr",        # For filtering transcript data
  "ArgumentCheck" # For validating function arguments
  "testthat",     # For running tests
  "roxygen2",     # For generating documentation
  "devtools",     # For installing chimeraviz
  "knitr"         # For creating vignettes
)
install.packages(cranPackages)
```

Then install these Bioconductor dependencies:

```
# Bioconductor
biocPackages <- c(
  "org.Hs.eg.db",      # For mapping gene names to ensembl ids
  "Rsamtools",         # For loading data from .bam files
  "GenomeInfoDb",      # For translating chromosome names (seqnames)
  "GenomicAlignments", # For loading data from .bam files
  "Biostrings",        # For representing sequences
  "GenomicRanges",     # For representing ranges
  "IRanges",           # For representing intervals
  "S4Vectors",         # For vectors
  "AnnotationDbi",     # For loading annotation files
  "ensembldb",         # For loading transcript data
  "AnnotationFilter",  # For filtering transcript data
  "graph",             # For plotting
  "Gviz",              # For plotting
  "Rgraphviz",         # For plotting
  "BiocStyle"          # For creating vignettes
)
source("https://bioconductor.org/biocLite.R")
biocLite(biocPackages)
```

And finally you can install `chimeraviz` with this command:

```
devtools::install_github(
  "stianlagstad/chimeraviz",
  build_vignettes = TRUE)
```

Please [create an issue on Github](https://github.com/stianlagstad/chimeraviz/issues) if you have any problems at all.

# Code examples

Run `browseVignettes("chimeraviz")` to see package vignette, which shows most of what `chimeraviz` can do. Also see the reference manual.

# Tests

Tests are written with [testthat](https://cran.r-project.org/web/packages/testthat/index.html) and are located in `tests/testthat`. They can be run with `devtools::test()` if you have cloned this repository, _i.e._ not installed the package with `devtools::install_github()` but have used `git clone git@github.com:stianlagstad/chimeraviz.git chimeraviz`.

# Credits

This package was developed by Stian Lågstad for his master thesis: Visualizing chimeric RNA. The work was supervised by [Rolf Skotheim](http://ous-research.no/skotheim/) and [Ole Christian Lingjærde](http://www.mn.uio.no/ifi/personer/vit/ole/).

# Licence

[Artistic Licence 2.0](https://opensource.org/licenses/Artistic-2.0).
