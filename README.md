# Discovery of cell subpopulations using sigle cell sequencing data and ReSET

# What is ReSET (Robust Subpopulation dEcision Trees)?

ReSET is a tool for discovering robust populations derived from single cell RNA-Seq experiments.

# Overview

Single-cell RNA sequencing enables unbiased analysis of expression patterns, but researchers don’t have the tools for appropriate decision making during the analysis. Our general aim is to introduce unbiased data-driven strategies to identify the appropriate number of robust subpopulations, their discriminatory defining markers and the relationships between populations.

# Workflow diagram


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
install_github("kbroman/broman")
```

binaries for ReSET are available through GNU Guix, the
functional package manager for reproducible, user-controlled software
management.  Install the complete pipeline bundle with the following
command:

```sh
guix package -i ReSET
```
# Example Data Sources
[10X Genomics, 4k Pan T Cells from a Healthy Donor](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/t_4k)

# Dependencies

python = 2.7
R >= 3.2

# Contributors (Alphabetical Order)
* [Assaf Magen]() (Team Lead)
* [Jonathan Badger]()
* [Billy Kim]()
* [Claire Malley]()
* [Chris Rhodes]()
* [Amulya Shastry]()
* [Mamie Wang]()


# License
