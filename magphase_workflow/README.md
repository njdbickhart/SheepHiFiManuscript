# MAGPhase workflow
---

This is a collection of useful scripts and a snakemake workflow designed to automate MAGPhase haplotyping. 

## MAGPhase 

This is a snakemake workflow designed to run the mag_phaser.py script that is part of the [cDNA_cupcake](https://github.com/Magdoll/cDNA_Cupcake) repository. The workflow will align HiFi reads to your assembly, phase the reads within target single copy gene regions, estimate the number of haplotypes and then generate tabular output summaries for each MAG. 

---

### Requirements

* Python 3.6+
* miniConda 3.6+

### Installation

To ensure reproducible output, we have provided a "yaml" file that can be used to create a virtual environment suitable to run MAGPhase. You will still need to install cDNA_cupcake as part of that environment before using the workflow. Here is an example of how to accomplish this using the [miniConda](https://docs.conda.io/en/latest/miniconda.html) package manager:

```bash
# This creates the environment
conda env create -f magphase.yaml

# Now you must activate the environment
conda activate magphase

# Finally, download and install cDNA_cupcake 
git clone https://github.com/Magdoll/cDNA_Cupcake.git
cd cDNA_Cupcake
python setup.py build
python setup.py install

# You should be all set! Remember to activate the environment before running this workflow in the future!
```
---
### Input

MAGPhase requires the following files as input to the workflow [input methods in brackets].

#### General requirements for the workflow

* The entire metagenome assembly fasta file [default.json]
* Your high-quality MAGs [in fasta format, in a "mags" folder]
* Single copy gene loci [in BED format, in a "beds" folder]
* Your HiFi reads [default.json]

Similar to most snakemake workflows, MAGPhase recognizes files in the current active directory using regular expression matching. Specifically, MAGPhase looks for files in the following two folders in your current directory.

#### Required formatted files

* mags/*.fa|sta
* beds/*.bed
* default.json

With the exception of the "default.json" file, these files must have a common prefix (in lieu of the asterisk), with each "mags" file having a corresponding "beds" file and vice versa. Here is an example of paired output suitable for running in the pipeline:

```bash
mags/bin356.fa
beds/bin356.bed
```

If your files are not suitably paired, MAGPhase will complain about missing input files!

Additionally, you must create a "default.json" file that points to other files relevant to the workflow. We have included an example file with this repository, so the best way to configure your workflow is to copy this template to your active directory and modify it in situ like so:

```bash
cp ~/MAGPhase/default.json ./
vim default.json
```

Replace text with settings specific to your environment as needed.

---

### Usage

After preparing all input and creating your virtual environment, you are ready to run MAGPhase! To start the workflow, navigate to the directory that contains the "mags" and "beds" folders and run the following:

```bash
snakemake -s SheepHiFiManuscript/magphase_workflow/MAGPhase --jobs 100
```

If you are running tasks on a cluster (preferred), then you can specify resource requirements for your cluster environment using the provided **cluster.json** file in this directory. Here is an example of how to specify cluster parameters in a Slurm environment:

```bash
snakemake -s SheepHiFiManuscript/magphase_workflow/MAGPhase --cluster-config SheepHiFiManuscript/magphase_workflow/cluster.json --cluster "sbatch -N 1 -n {cluster.threads} --mem {cluster.mem} -o {cluster.stdout} -p priority" --jobs 100
```

Note that you can override the submission string parameters with static values, or you can copy the cluster.json file and modify the contents so that they are more suitable to your environment.

---

### Output

MAGPhase will generate four folders in your working directory and each will contain specific output files from the entire workflow:

| Folder name | Sub files/folders | Description |
| :--- | :--- | :--- |
| **logs** | - | Runtime logging information for MAGPhase. Check these for status reports on each rule |
| **mapping** | - | This folder will contain a bam file with HiFi read mappings to your assembly fasta file. |
| **phase** | - | These folders contain the results of the main MAGPhase workflow. The "MAG" prefix will be replaced by the specific MAG name that was analyzed in this task. |
| ^ | MAG.NO_SNPS_FOUND	| Only present when no SNPs were identified! Indicates that this MAG has no identifiable SNP variants suitable for haplotyping |
| ^ | MAG.human_readable_by_hap.txt | Each identified SNP haplotype is listed on each line, along with average variant read counts (count) that support the entire haplotype |
| ^ | MAG.human_readable_by_pos.txt | Provides exact variant allele counts (count) and zero-based genomic coordinates (contig, pos) for each variant site within the haplotype (indicated by the "varIdx" column) in descending order. |
| ^ | MAG.human_readable_by_read.txt | Lists all reads and their correspondence to identified haplotypes. Very useful for color-coding reads in viewers such as [IGV](https://software.broadinstitute.org/software/igv/) |
| **consolidate** | - | These folders contain merged results from the three specific output files produced during the **phase** step of the workflow. The files are useful for plotting and tabulating general statistics from your assembly. |
| ^ | consolidated.short | A file that lists the maximum number of observed haplotypes (column 2) for each MAG(column 1). A count of "0" indicates that No suitable SNP haplotypes were identified and that the MAG appears to be lineage-resolved |
| ^ | consolidated.long | This file provides average read depth estimates for each observable haplotype string within user supplied BED coordinate regions. If the MAG did not have identified haplotypes in a region, it is not listed |

Column descriptions for **consolidated.long**:

1. The MAG name
2. The SNP haplotype string (only contains variant SNPs)
3. The length of the SNP haplotype
4. The average number of reads supporting the SNP haplotype.
5. The contig name for the region
6. The start of the region
7. The end of the region


If you would like to visualize the results of your run, we have provided a script (**plotMagPhaseOutput.py**, see below) to generate read depth plots of each SNP haplotype in your MAG. 

---

## Scripts

Two scripts are provided to assess the quality of your MAG data and to visualize the output of MAGPhase. 

---

### getBAMMapQ0Ratios.py

This script performs a window analysis on a BAM file to generate a ratio of MapQ0 reads (ambiguously mapped reads) to total read depth. This is useful to assess the likelihood of mapping shorter reads to your assembled MAGs.

#### Usage:

```bash
usage: getBAMMapQ0Ratios.py [-h] -b BAM -o OUTPUT [-w WINDOWLENGTH]

Window analysis to calculate ratio of mapq0 reads in bam files

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     Input bam file. Must be coord sorted and indexed.
  -o OUTPUT, --output OUTPUT
                        Output file name. Output format is bed format
  -w WINDOWLENGTH, --windowlength WINDOWLENGTH
                        Window length for contig windows [5000]
```

Output bed files have the following columns:

1. Contig Name
2. Start bp
3. End bp
4. Boolean flag if the window contains any MapQ0 reads
5. Ratio of MapQ0 reads to total reads
6. Total read count
7. MapQ0 read count

---

### plotMagPhaseOutput.py

This script uses the python seaborn library to plot a linear graph of read depth for MAGs characterized by the MAGPhase workflow. 

#### Usage:

```bash
usage: plotMagPhaseOutput.py [-h] -f FAI -o OUTPUT -b BAM -u HUMAN
                             [-i BINSIZE] [-e BREAKS]

A tool to plot bin and contig level read depth differences in strain
assignment

optional arguments:
  -h, --help            show this help message and exit
  -f FAI, --fai FAI     Input reference fasta index file for the bin
  -o OUTPUT, --output OUTPUT
                        Output file basename. Output files are {output}.wins
                        and {output}.pdf
  -b BAM, --bam BAM     Input CCS read depth bam file
  -u HUMAN, --human HUMAN
                        Input human-readable position variant call file
  -i BINSIZE, --binsize BINSIZE
                        Bin size in bases [5000 bp]
  -e BREAKS, --breaks BREAKS
                        Draw windows in the following bed regions [Optional]
```

#### Arguments to the script

| Argument | Optional? | Description |
| :---   | :--- | :---- |
| -f | No | A fasta index file for the input MAG (can be generated with the *samtools faidx* command |
| -b | No | A sorted and indexed BAM file containing alignments of HiFi reads to the assembly. Does not need to be filtered to contigs specifically in the MAG, but must contain alignments to all contigs in the MAG |
| -u | No | The *MAG.human_readable_by_pos.txt* text file generated in the *phase* folder of the MAGPhase workflow for this specific MAG. |
| -o | No | The output file prefix. Two files will be generated with this prefix, and both are described below. |
| -i | Yes | Designate the size of the window in base pairs (integer). Larger window sizes smooth out random fluctuations of coverage, but be careful not to extend your window size larger than your smallest single copy gene region used for the MAG in the workflow! |
| -e | Yes | To annotate the plot, you can add a bed file (only the first three fields are required) to draw boundaries of single copy gene regions, or regions of other importance. |

#### Output files:

This script generates one text file (*output.wins*) and one image file (*output.pdf*) when it successfully concludes. 

The *.wins* file is a tab-delimited table with a header. The "hap" column indicates either the reference depth of coverage (for "0") or the alternate haplotype rank (sorted by read depth). The "count" column gives the depth of read coverage for the reference or haplotype in the window. 

the *.pdf* file is a seaborn plot of the total MAG length, with vertical black lines indicating the boundaries of contigs in the MAG. The reference depth of coverage is indicated by a blue line, whereas alternative haplotype depths of coverage are represented by alternating colors. 

**Note** that the "reference" color will drop in regions of alternative haplotype discovery. This is expected, as shared regions of the genome among lineages will have higher depths of coverage when not separated out into separate, lineage-resolved MAGs. 