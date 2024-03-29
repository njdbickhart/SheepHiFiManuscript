# SheepHiFiManuscript

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5120910.svg)](https://doi.org/10.5281/zenodo.5120910)

A collection of scripts and workflows to analyze HiFi assembled metagenomes. Each folder contains workflows and scripts designed to interrogate metagenome assembly data as demonstrated in our [preprint publication](https://www.biorxiv.org/content/10.1101/2021.05.04.442591v1.abstract). 

## Installation

To clone this repository and the associated MAGPhase program, run the following command:

```bash
git clone --recursive https://github.com/njdbickhart/SheepHiFiManuscript.git
```

Each of the following folders has separate installation instructions and runtime requirements, so please consult the provided README files for more details. We tried to limit the use of external libraries and programs where possible, so these requirements shouldn't be too onerous! 

## Downsampling

This folder contains scripts and a workflow to replicate our progressive downsampling analysis of HiFi sequence data. 

We have provided a script to generate downsampled read subsets from HiFi data, and an example snakemake workflow to replicate the analysis conducted in the manuscript. 

For more details, please consult the README.md file in that directory.

## MAGPhase_workflow

This folder contains the packaged workflow for the MAGPhase algorithm in the form of a snakemake pipeline. Additionally, we provide a script to visualize the output of the workflow as a separate python script. 

We have provided detailed instructions on how to install and run this workflow in the README.md file in the folder.

## Viral_association

This folder contains an updated, reposted pipeline script for our viral association analysis. We previously posted this pipeline as part of our publication on the [cattle rumen metagenome in Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1760-x) and ask that you cite that manuscript if you use this pipeline. New features added to this version of the pipeline include the automated generation of network plots to visualize the association of viral contigs with candidate host bacteria and archaea.

More details are provided in the README.md file in this directory.