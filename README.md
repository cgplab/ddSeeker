# ddSeeker
A tool for processing Bio-Rad ddSEQ single cell RNA-seq data.

# Description
`ddSeeker` identifies cellular and molecular identifiers from single cell RNA sequencing experiments.

Users must provide either one unmapped BAM file
([uBAM](https://gatkforums.broadinstitute.org/gatk/discussion/11008/ubam-unmapped-bam-format))
or two fasta files from a paired-end read experiment, and a text file with a
set of barcode blocks, one per line (if you don't have such file
`ddSeeker_barcodes.py` can generate one).

The program returns a new unmapped BAM file with tagged reads. Default tags are
**XB** and **XU** for cell and molecular identifiers respectively, and **XE**
for different errors related to barcode identification.

##### Additional options

  - Work with multiple cores using `-c/--cores`.
  - Generate 2 summary files reporting statistics for cell barcodes and errors
    of tag retrival using `-s/--summary-prefix`.
  - Output to standard output (allowing direct feeding to other tools for
    filtering, sorting etc.).
  - Manually set [tags](https://genome.sph.umich.edu/wiki/SAM#What_are_TAGs.3F)
    with `--tag-bc`, `--tag-umi` and `--tag-error` options.

## Examples
###### Generate barcodes file<sup>1</sup>

    ddSeeker_barcodes.py -i sampleA_R1.fastq.gz -c 20

###### Run ddSeeker with one unmapped bam

    ddSeeker.py -i sampleA_unmapped.bam -o sampleA_tagged.bam -c 20

###### Run ddSeeker with fastq files

    ddSeeker.py -i sampleA_R1.fastq.gz sampleA_R2.fastq.gz -o sampleA_tagged.bam -c 20

###### Run ddSeeker, print to stdout (use "-") and pipe to samtools for queryname sorting

    ddSeeker.py -c 20 -i sampleA_R* -o - | samtools sort -n -o sampleA_tagged_qsorted.bam

###### Run ddSeeker and produce summary files

    ddSeeker.py -i sampleA_unmapped.bam -o sampleA_tagged.bam -s summary_folder/sampleA -c 20

###### Run full ddSeeker-"drop-seq tools" pipeline
    ddSeeker_alignment.sh [options] sample1_unmapped.bam

## Dependencies
- [Python](https://www.python.org/downloads/release/python-365/) (>= 3.5)
- [Biopython](http://biopython.org/wiki/Download) (>= 1.71)
- [pysam](https://pysam.readthedocs.io/en/latest/index.html) (>= 0.14)
<!-- - [Drop-seq_tools](http://mccarrolllab.com/dropseq/) (>= 1.13) -->
<!-- - [STAR](https://github.com/alexdobin/STAR) (>= 2.5) -->

### Installing dependencies
We suggest to use `pip` to install python packages which should be installed
by default with Python > 3.5.

######

    pip install biopython
    pip install pysam


###### Footnotes
<sup>1</sup>Barcodes are extracted from the most represented strings at fixed
positions. Ideally, all the bam/fastq files from a single cell RNA-seq
experiment should be processed. Practically, even a single file would be
sufficient.
