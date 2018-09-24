FROM rocker/verse:3.5.1
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('limma')"
