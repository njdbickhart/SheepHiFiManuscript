# Downsampling analysis
---

This is a collection of scripts designed to downsample HiFi reads for assembly analysis. The scripts and associated workflow are provided in the hopes that they will be useful, but are otherwise provided without guarantee of extensive troubleshooting.

Scripts have only been tested on local hardware and are designed to work on distributed computing systems using the miniconda package manager to generate virtual environments. 

## downsampleFastaForAssembly.pl

This script is designed to recursively downsample a set of HiFi read fasta files using a progressive sampling strategy. Briefly, Individual reads are randomly sampled from each prior strata in descending order to generate individual downsampled fasta files. In other words, reads in the "10%" strata must all be present in the "90%" strata as part of this strategy. 

#### Requirements

This script only requires base Perl v5.8+ without any additional modules.

#### Usage:

```bash
perl downsampleFastaForAssembly.pl <fasta file 1> ... <fasta file n>
```

#### Output:

As written, the script will generate nine subsampled fasta files, each named "downsample.X.fasta," with the "X" denoting the percentage of reads sampled from the input. The fasta headers and sequence will be otherwise unaltered, apart from the newline formatting of the actual sequence being removed in the downsampled files.

We recommend assembling downsampled reads for further analysis using the "--pacbio-hifi" input argument in [metaFlye](https://github.com/fenderglass/Flye).

## Downsample

This is a snakemake workflow that was used to automate the binning and analysis of downsampled reads. We encourage other groups to co-opt sections of this code for use on their projects, but many configuration settings will need to be modified within the file itself to make it work for other environments.

If you wish to run this workflow on your system, please edit the "config" dictionary on line 11 of the snakefile to point to the files, executables and settings that comply to your dataset.

#### Requirements

You must have the following installed at a minimum:

* Python 3.6+
* Snakemake 5.8+
* miniconda 3.6+

#### Usage:

Make sure that you have properly configured the workflow to your environment first! To run the workflow, type the following:

```bash
snakemake -s Downsample --use-conda --jobs 50 --verbose --latency-wait 40
```


 