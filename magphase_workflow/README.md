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

MAGPhase requires the following files as input to the workflow [input methods in brackets]:

* The entire metagenome assembly fasta file [config.json]
* Your high-quality MAGs [in fasta format, in a "mags" folder]
* Single copy gene loci [in BED format, in a "beds" folder]
* Your HiFi reads [config.json]

Similar to most snakemake workflows, MAGPhase recognizes files in the current active directory using regular expression matching. Specifically, MAGPhase looks for files in the following two folders in your current directory:

* mags/*.fa|sta
* beds/*.bed
* config.json

With the exception of the "config.json" file, these files must have a common prefix (in lieu of the asterisk), with each "mags" file having a corresponding "beds" file and vice versa. Here is an example of paired output suitable for running in the pipeline:

```bash
mags/bin356.fa
beds/bin356.bed
```

If your files are not suitably paired, MAGPhase will complain about missing input files!

Additionally, you must create a "config.json" file that points to other files relevant to the workflow. We have included an example file with this repository, so the best way to configure your workflow is to copy this template to your active directory and modify it in situ like so:

```bash
cp ~/MAGPhase/config.json ./
vim config.json
```

Replace text with settings specific to your environment as needed.

---

### Usage



## Scripts

Two scripts are provided to 