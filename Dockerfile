FROM rocker/verse:latest
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('limma')"
