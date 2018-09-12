# Discovery of cell subpopulations using single cell sequencing data and ReSET

# What is ReSET (Robust Subpopulation dEcision Trees)?

ReSET is a tool that leverages the knowledge of cell cluster heirarchy to discover robust populations derived from single cell RNA-Seq experiments. ReSET creates a heirarchical model that fits two independent single cell RNA-Seq datasets.

# Overview

Single-cell RNA sequencing enables unbiased analysis of expression patterns, but researchers don’t have the tools for appropriate decision making during the analysis. Our general aim is to introduce unbiased data-driven strategies to identify the appropriate number of robust subpopulations, their discriminatory defining markers and the relationships between populations.

# Workflow diagram
![Alt text](https://github.com/NCBI-Hackathons/robustSingleCell/blob/master/ReSET%20pipeline.png?raw=true)

# Getting started

To run ReSET on your experimental data, describe your samples in a CSV
file `sample_sheet.csv`, provide a `settings.yaml` to override the
defaults defaults, and select the pipeline.

To generate a settings file template for any pipeline:
```sh
ReSET [pipeline] --init=settings
```

To generate a sample sheet template for any pipeline:
```sh
ReSET [pipeline] --init=sample-sheet
```

Here's a simple example to run the RNAseq pipeline:\
```sh
ReSET rnaseq my-sample-sheet.csv --settings my-settings.yaml
```

To see all available options run `ReSET --help`.

# Install
A pre-built package is available in this repository.

First, you need to install the [devtools](https://github.com/hadley/devtools) package which is available from CRAN. Invoke R and then type

```R
install.packages("devtools")
```

Load the devtools package.
```R
library(devtools)
```

Install ReSET
```R
install_github("NCBI-Hackathons/robustSingleCell")
```

# Example Data Sources
[10X Genomics, 4k Pan T Cells from a Healthy Donor](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/t_4k)

[10X Genomics, 3k Pan T Cells from a Healthy Donor](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/t_3k)

# Dependencies

python = 2.7

R >= 3.5

# Contributors (Alphabetical Order)
* [Assaf Magen](https://github.com/asmagen) (Team Lead)
* [Jonathan Badger](https://github.com/jhbadger)
* [Billy Kim](https://github.com/bkim62)
* [Claire Malley](https://github.com/cemalley)
* [Chris Rhodes](https://github.com/ctrhodes)
* [Amulya Shastry](https://github.com/amulyashastry)
* [Mamie Wang](https://github.com/Mamie)


# License
