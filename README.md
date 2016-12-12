# chimeraviz

chimeraviz is an R package that automates the creation of chimeric RNA visualizations.

# Installing chimeraviz

`chimeraviz` will be submitted to Bioconductor, but until then the package can be installed manually. First install CRAN dependencies:

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
  "plyr"          # For converting fusion data to data.frame format
)
install.packages(cranPackages)
```

Then install Bioconductor dependencies:

```
# Bioconductor
biocPackages <- c(
  "org.Hs.eg.db",      # To map gene names to ensembl ids
  "Rsamtools",         # To load data from .bam files
  "GenomeInfoDb",      # To translate chromosome names (seqnames)
  "GenomicAlignments", # To load data from .bam files
  "Biostrings",        # To represent sequences
  "GenomicRanges",     # To represent ranges
  "IRanges",           # To represent intervals
  "S4Vectors",         # To get vectors
  "AnnotationDbi",     # To load annotation files
  "GenomicFeatures",   # To load data from annotation databases
  "graph",             # To plot graphs
  "Gviz",              # To plot various plots
  "Rgraphviz",         # To plot graphs
  "BiocStyle"          # To create vignettes
)
source("https://bioconductor.org/biocLite.R")
biocLite(biocPackages)
```

You can then install `chimeraviz` with this command:

```
devtools::install_github(
  "stianlagstad/chimeraviz",
  build_vignettes = TRUE)
```

Please create an issue if you have any problems at all.

# Code examples

Please run `browseVignettes("chimeraviz")` code examples in the package vignette.

# Tests

Tests are written with [testthat](https://cran.r-project.org/web/packages/testthat/index.html) and are located in `tests/testthat`. They can be run with `devtools::test()` if you have cloned this repository, _i.e._ not installed the package with `devtools::install_github()` but have used `git clone git@github.com:stianlagstad/chimeraviz.git chimeraviz`.

# Credits

This package was developed by Stian Lågstad for his master thesis: Visualizing chimeric RNA. The work was supervised by [Rolf Skotheim](http://ous-research.no/skotheim/) and [Ole Christian Lingjærde](http://www.mn.uio.no/ifi/personer/vit/ole/).

# Licence

[Artistic Licence 2.0](https://opensource.org/licenses/Artistic-2.0).
