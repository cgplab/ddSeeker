## ddSeeker
A tool for processing single cell RNA-seq ddSEQ data.

### Details
ddSeeker identifies cellular and molecular identifiers from single cell RNA sequencing experiments.

### Input
- Unmapped BAM file
  [uBAM](https://gatkforums.broadinstitute.org/gatk/discussion/11008/ubam-unmapped-bam-format)
  from a paired-end read experiment: FASTQ files can be converted using
  [picard FastqToSam](https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam)
- Text file with a set of barcode blocks, one per line (if not available, run `ddSeeker_barcodes.py` to generate one).

### Output
- Tagged [tagged](https://genome.sph.umich.edu/wiki/SAM#What_are_TAGs.3F) uBAM
  file with either the correct identifiers ("XB", "XU"), or an error ("XE") to specify
  why the identification has failed.
- Optionally, use `-s/--summary` to generate a summary of unique barcodes
  and error count.

### Example
    java -jar picard.jar FastqToSam F1=read1.fastq.gz F2=read2.fastq.gz O=output.bam
    ddSeeker\_barcodes.py \*.bam
    ddSeeker.py sample1.bam sample1\_tagged.bam -s summary/sample1 -n 10

### Dependencies
- [Python](https://www.python.org/downloads/release/python-365/) (>= 3.5)
- [Biopython](http://biopython.org/wiki/Download) (>= 1.71)
- [pysam](https://pysam.readthedocs.io/en/latest/index.html) (>= 0.14)
- [Drop-seq\_tools](http://mccarrolllab.com/dropseq/) (>= 1.13)
