# ddSeeker
A tool for processing Bio-Rad ddSEQ single cell RNA-seq data.

### Description
**ddSeeker** identifies cellular and molecular identifiers from single cell RNA sequencing experiments.

**Input**: R1 and R2 FASTQ files from a paired-end single cell sequencing experiment.

**Output**: one [unmapped BAM](https://gatkforums.broadinstitute.org/gatk/discussion/11008/ubam-unmapped-bam-format)
file containing reads tagged with cell barcodes and unique molecular identifiers (UMI).
Default [tags](https://genome.sph.umich.edu/wiki/SAM#What_are_TAGs.3F) are
**XB** and **XU** for cell barcodes and UMI, and **XE** for errors related to
the barcode identification.
Users can manually set different tags (see #additional-options).

#### Errors in barcode identification

    - LX = both linkers not aligned correctly
    - L1 = linker 1 not aligned correctly
    - L2 = linker 2 not aligned correctly
    - I  = indel in BC2
    - D  = deletion in Phase Block or BC1
    - J  = indel in BC3 or ACG trinucleotide
    - K  = indel in UMI or GAC trinucleotide
    - B  = one BC with more than 1 mismatch

#### Additional options

  - Increment number of CPU units (faster analysis) with `-c/--cores`.
  - Manually set tags with `--tag-bc`, `--tag-umi` and `--tag-error`.
  - Print uncompressed SAM file to standard output (allowing direct feeding to other tools for
    filtering, sorting etc.) with `-o/--output -` (note the `-` sign).
  - Generate two csv files reporting the number of reads per cell and the distribution
  of error tags specifying the path with `-s/--summary-prefix`.
  - Create plots from the csv summary files using `make_graphs.R` (see ).


### Install ddSeeker
Clone the repository and add the folder to your PATH variable

    git clone https://github.com/cgplab/ddSeeker.git
    export PATH=<path_to_ddSeeker>:$PATH

### Dependencies
- [Python](https://www.python.org/downloads) (>= 3.5)
- [Biopython](http://biopython.org) (>= 1.71)
- [pysam](https://pysam.readthedocs.io) (>= 0.14)

We suggest to install python packages using [pip](https://pip.pypa.io/en/stable/installing/)
which should be already installed if you are using Python3 >= 3.4.

    pip install -U pip  # upgrade pip
    pip install biopython
    pip install pysam

## Usage examples

##### ddSeeker with 20 cores

    ddSeeker.py --input sampleA_R1.fastq.gz sampleA_R2.fastq.gz --output sampleA_tagged.bam --cores 20

##### Print to stdout and pipe to samtools for queryname sorting

    ddSeeker.py -i sampleA_R* -c 20 -o - | samtools sort -no sampleA_tagged_qsorted.bam

##### Generate summary files and make graphs
Requires [R >=3.4](https://www.r-project.org/) and the [tidyverse](https://www.tidyverse.org/) package.
Three plots are generated: dot plot of error distribution, absolute count of reads per
cell, and cumulative distribution of reads per cell. The latter two report by default
the whole set of barcodes in the csv file. To limit the report to a lower
number, specify it from the command line.

    mkdir summary_folder
    ddSeeker.py -i sampleA_R* -c 20 -o sampleA_tagged.bam -s summary_folder/sampleA
    make_plot.R summary_folder/sampleA 2000

## Integrating single cell analysis pipelines
Several pipelines have been developed to perform single cell analysis.
Below we describe the main steps required to integrate our tool with
[Drop-seq tools](http://mccarrolllab.com/dropseq/),
[scPipe](https://github.com/LuyiTian/scPipe) and
[dropEst](https://github.com/hms-dbmi/dropEst).


##### Drop-seq tools
Since Drop-seq tools was our choice for our analyses, we provide a ready-to-use bash
script.  Simply run

    ddSeeker_dropSeq_tools.sh [options] sampleA_R1.fastq.gz sampleA_R2.fastq.gz

to produce aligned tagged reads in BAM format. 
Table of Counts can be obtained using the `DigitalExpression` tool included in Drop-seq tools.
To do that, run it setiing `CELL_BARCODE_TAG=XB` and `MOLECULAR_BARCODE_TAG=XU`.

##### scPipe
scPipe requires one FASTQ file with cell barcodes and UMIs stored in the header
of each read record. To change the output of **ddSeeker** use the option
`--pipeline scPipe`.

    ddSeeker.py -i sampleA_R* -o sampleA_tagged.fastq.gz -c 20 --pipeline scPipe

In addition, set `bc_len=18` and `UMI_len=8` with the `sc_exon_mapping()` function.

##### dropEst
dropEst can work with tagged BAM files. Simply make the BamTags match with the
ddSeeker tags specifying them in the config.xml file

    <BamTags>
        <cb>XB</cb>
        <umi>XU</umi>
    </BamTags>

## Citation
Romagnoli et al., **ddSeeker: a tool for processing Bio-Rad ddSEQ single cell RNA-seq data**, submitted
