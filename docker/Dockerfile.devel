from bioconductor/devel_core2

RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("data.table")'
RUN R -e 'install.packages("roxygen2")'
RUN R -e 'install.packages("RCircos")'
RUN R -e 'install.packages("DT")'
RUN R -e 'install.packages("dplyr")'
RUN R -e 'install.packages("ArgumentCheck")'
RUN R -e 'install.packages("testthat")'
RUN R -e 'install.packages("lintr")'
RUN R -e 'install.packages("gtools")'

RUN R -e 'source("https://bioconductor.org/biocLite.R");biocLite("org.Hs.eg.db")'
RUN R -e 'source("https://bioconductor.org/biocLite.R");biocLite("org.Mm.eg.db")'

