# Viral Analysis Scripts
![DOI](https://zenodo.org/badge/159530113.svg)
---
This workflow was released as part of our first [Rumen Metagenome assembly manuscript](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1760-x) but is reproduced here for convenience.

## Viral (or mobile DNA) link analysis

For ease of use, I have encapsulated the workflow of our rumen metagenome analysis manuscript into a single script (**viralAssociationPipeline.py**) that automates the entire process! The wrapper is flexible, and will produce association tables for individual long-read datasets, Hi-C reads or merge the two datasets together. Here is a brief description of the pipeline and the requirements to run it:

### Usage statement
---

```
python3 viralAssociationPipeline.py  --help
usage: viralOverhangPort.py [-h] [-l LONG_READ] -a ASSEMBLY [-g VIRAL_CONTIGS]
                            -b BLOB_TOOLS [-i HIC_LINKS] -v VIRUSES
                            [-c LINK_THRESH] [-e OVERHANG] [-n] -o OUTPUT
                            [-m MINIMAP] [-s SAMTOOLS]

Process long-read alignments and Hi-C links to identify likely viral-host
associations in metagenomic assembly data

optional arguments:
  -h, --help            show this help message and exit
  -l LONG_READ, --long_read LONG_READ
                        Input long-read fastq file
  -a ASSEMBLY, --assembly ASSEMBLY
                        Fasta of the assembled contigs (in a single file) for
                        analysis
  -g VIRAL_CONTIGS, --viral_contigs VIRAL_CONTIGS
                        Fasta file of the separated contigs of viral origin (a
                        subset of the full assembly)
  -b BLOB_TOOLS, --blob_tools BLOB_TOOLS
                        Input blob tools or taxonomic data table
  -i HIC_LINKS, --hic_links HIC_LINKS
                        Input sam/bam file with alignments of Hi-C reads
  -v VIRUSES, --viruses VIRUSES
                        Tab delimited list of contigs containing viral
                        sequence and their lengths
  -c LINK_THRESH, --link_thresh LINK_THRESH
                        Filter for the number of stdevs of Hi-C links above
                        the average count to be used in viral Hi-C association
                        [2.5]
  -e OVERHANG, --overhang OVERHANG
                        Filter for long-read overhang off of viral contigs
                        [150]
  -n, --noplot          [optional flag] Disables plotting of data
  -o OUTPUT, --output OUTPUT
                        Output basename
  -m MINIMAP, --minimap MINIMAP
                        Path to the minimap executable
  -s SAMTOOLS, --samtools SAMTOOLS
                        Path to the samtools executable

```

### Requirements
---

* Python 3.6+
* Samtools 1.9+
* Minimap2
* NetworkX 

### Input datasets
---

Please note that you **MUST** input at least one Long Read (**-l, --long_read**) dataset or one Hi-C (**-i, --hic_links**) dataset to run the pipeline. All other arguments are required. 

| Name | Argument | Description |
| :--- | :--------| :-----------|
| Virus Names | **-v, --viruses** | A tab delimited file giving the viral contig name and the length of the contig. If you use *samtools faidx* on a separate fasta of your virus contigs you can pass the *.fai* file to this argument |
| Viral Contigs | **-g, --viral_contigs** | A separate fasta file of your viral contigs. *OPTIONAL* if not supplied, the program will try to generate this from the assembly fasta you entered |
| Assembly | **-a, --assembly** | A fasta file containing your assembled contigs (should also include your viral contigs!) |
| Blob Tools | **-b, --blob_tools** | This is a taxonomic file that can be generated via the [Blob Tools taxify command](https://blobtools.readme.io/docs/taxify). If you have an aversion to Blob Tools (why??) then there is an alternative [below](#blobtools) |
| Long Reads | **-l, --long_read** | This is a fasta file containing your long-reads (preferrably HiFi or low-error ONT) for alignment to the assembly. *OPTIONAL* If not included, this analysis will be skipped |
| Hi-C Alignments | **-i, --hic_links** | This is a SAM/BAM file containing the alignment of paired-end Hi-C reads to your assembly fasta. *OPTIONAL* If not included, this analysis will be skipped |
| Output Name | **-o, --output** | OK, this is not exactly an "input parameter" but it is required. Just provide a name and all of your output files will be prefixed by that same name! |

### Input parameters
---

I've tried to set default parameters for each of these settings, but you may need to tweak them to deal with aberrations in your sample. 

| Name | Argument | Default | Description |
| :--- | :--------| :-------| :-----------|
|Link Threshold | -c, --link_thresh | 2.5 (float) | The number of standard deviations above the mean to use for identifying significant Hi-C link counts in your **-i** dataset. Higher is more stringent. |
|Overhang | -e, --overhang | 150 (int) | The number of bases that your **-l** dataset reads must overlap viral contigs and host contigs. Higher is more stringent.|
|Minimap | -m, --minimap | "minimap2" (str) | Path to the minimap2 executable. |
|Samtools| -s, --samtools| "samtools" (str) | Path to the samtools executable. |
|Plotting| -n, --noplot | N/A | Disables investigative plotting using the python3 networkX module |

### Output files
---

Just to document progress and allow you to interrogate the data at each step, the script will produce the following output files (all tab-delimited; in order):

* **(OUTPUT).algn.viruses** : Contains the overhang coordinates for the **-l** dataset on a per-read basis. Columns:
	1. Read name
	2. Start
	3. End
	4. Viral Contig
* **(OUTPUT).lread.vir.graph**: Contains a di-graph table of viral-host links from the **-l** dataset on a per-read basis.
* **(OUTPUT).hiclinks.tab** : Consists of a count (column 3) of Hi-C links between viral (column 1) and candidate host (column 2) contigs.
* **(OUTPUT).final.tab** : What you came here for! The final table containing your association data! Columns:
	1. Viral contig
	2. Host contig
	3. Association {Read, HiC, or Both}
	4. Viral Genus
	5. Host Kingdom
	6. Host Genus
	7. Evidence (semi-colon delimited counts of reads and/or Hi-C links for this association)

<a name="blobtools"></a>
### Creating a taxonomic affiliation file (sans Blob tools)
---

For this pipeline, we are particularly interested in the putative taxonomic affiliation of each contig as assessed via the Blobtools [taxify](https://blobtools.readme.io/docs/taxify) workflow. If you have a strong aversion to Blobtools, or if it just doesn't fit in your workflow, you can simulate this file by generating a tab-delimited text file with the following information:

```
## Ignored comments are prefaced with double "hashes"
# Contig\tsuperkingdom.t\tgenus.t (single hash denotes the header. Must contain these three fields. Only the position of the "Contig" column must be in the first column of the file)
contig1\tBacteria\tPrevotella
contig2\tBacteria\tPseudomonas
```


