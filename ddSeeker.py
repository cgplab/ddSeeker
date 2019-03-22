#!/usr/bin/env python3

import logging
import pysam
import sys
import time
from argparse import ArgumentParser
from Bio import pairwise2
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from functools import partial
from gzip import open as gzopen
from itertools import islice
from multiprocessing import Pool
from numpy import cumsum
from pathlib import Path
from re import match as re_match

logging.basicConfig(level=logging.INFO, datefmt='%H:%M:%S',
                    format="[%(asctime)s] %(levelname)s - %(message)s")
_linkers = ["TAGCCATCGCATTGC", "TACCTCTGAGCTGAA"]
_local_aligner = partial(pairwise2.align.localxs, one_alignment_only=True)
_global_aligner = partial(pairwise2.align.globalxs, score_only=True, one_alignment_only=True)

_cell_count={}
_error_count = {}


def hamming_dist(s1, s2):
    """Return the Hamming distance between equal-length sequences

    >>> hamming_distance("ACG", "ACT")
    1
    """
    if (len(s1) != len(s2)):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

def fix_block(block):
    """Set barcode to most similar (up to 1 mismatch)

    >>> fix_block("TTTGGG")
    'TTTGGG'

    >>> fix_block("TTTGAG")
    'TTTGGG'

    >>> fix_block("TTTGAA")
    ''
    """
    for bc in _barcodes:
        score = _global_aligner(bc, block, -1, -1)
        if len(block) >= 5 and score >= 5:
            return(bc)
    else:
        return(None)

def get_tags(sequence):
    """Extract barcodes from R1 and return a tuple of SAM-format TAGs.
    Default tag values:
    - XC = cell barcode
    - XM = molecular barcode (UMI)
    - XE = error

    Errors:
    LX = both linkers not aligned correctly
    L1 = linker 1 not aligned correctly
    L2 = linker 2 not aligned correctly
    I = indel in BC2
    D = deletion in Phase Block or BC1
    J = indel in BC3 or ACG trinucleotide
    K = indel in UMI or GAC trinucleotide
    B = one BC with more than 1 mismatch"""

    sequence = sequence.upper()
    starts = []
    k = []
    for linker in _linkers:  # align the two linkers
        alignment = _local_aligner(sequence, linker, -2, -1)[0]
        seqA, seqB, score, begin, end = alignment
        length = end - begin
        if length == 15 and score >= 14:
            # 0-1 mismatch
            starts.append(begin)
            k.append(0)
        elif length == 14 and score == 14:
            # 1 mismatch at starting position
            starts.append(begin - 1)
            k.append(0)
        elif "-" in seqA[begin:end] and length == 15 and score == 12:
            # 1 deletion
            starts.append(begin)
            k.append(-1)
        elif "-" in seqB and length == 16 and score == 13:
            # 1 insertion
            starts.append(begin)
            k.append(1)
        else:
            starts.append(None)
            k.append(None)

    if not starts[0] and not starts[1]:
        return(dict([(_tag_error, "LX")])) # no linker aligned
    elif not starts[0]:
        return(dict([(_tag_error, "L1")])) # linker 1 not aligned
    elif not starts[1]:
        return(dict([(_tag_error, "L2")])) # linker 2 not aligned

    if starts[1]-starts[0] == 21+k[0]:
        bc2 = sequence[starts[1]-6: starts[1]]
    elif starts[1]-starts[0] == 20+k[0]: # 1 deletion in bc2
        bc2 = sequence[starts[1]-5: starts[1]]
    elif starts[1]-starts[0] == 22+k[0]: # 1 insertion in bc2
        bc2 = sequence[starts[1]-7: starts[1]]
    else:
        return(dict([(_tag_error, "I")]))

    if starts[0] < 5:
        return(dict([(_tag_error, "D")]))
    elif starts[0] == 5:
        bc1 = sequence[: starts[0]]
    else:
        bc1 = sequence[starts[0]-6: starts[0]]

    acg = sequence[starts[1]+21+k[1]: starts[1]+24+k[1]]
    i = 0
    thr = 0
    while i < len(acg) and thr <= 1:
        if acg[i] != "ACG"[i]:
            thr += 1
        i += 1
    else:
        if thr > 1:
            return(dict([(_tag_error, "J")]))

    try:
        dist_acg = hamming_dist(acg, "ACG")
    except ValueError:
        dist_acg = float("inf")
    if dist_acg > 1:
        return(dict([(_tag_error, "J")]))

    gac = sequence[starts[1]+32+k[1]: starts[1]+35+k[1]]
    try:
        dist_gac = hamming_dist(gac, "GAC")
    except ValueError:
        dist_gac = float("inf")
    if dist_gac > 1:
        return(dict([(_tag_error, "K")]))

    bc3 = sequence[starts[1]+15+k[1]: starts[1]+21+k[1]]

    barcode = []
    for block in (bc1, bc2, bc3):
        fixed = fix_block(block)
        if fixed:
            barcode.append(fixed)
        else:
            return(dict([(_tag_error, "B")]))

    umi = sequence[starts[1]+24+k[1]: starts[1]+32+k[1]]

    return(dict([(_tag_bc, "".join(barcode)), (_tag_umi, umi)]))

def compute_summary(tags):
    # summary statistics
    if tags.get(_tag_error):
        _error_count[tags[_tag_error]] = _error_count.get(tags[_tag_error], 0) + 1
    elif tags.get(_tag_bc):
        _cell_count[tags[_tag_bc]] = _cell_count.get(tags[_tag_bc], 0) + 1
        _error_count["PASS"] = _error_count.get("PASS", 0) + 1

def write_summary(summary_path):
    file_name = summary_path + ".errors.csv"
    ordered_tags = ["LX", "L1", "L2", "I", "D", "J", "K", "B", "PASS"]
    out = open(file_name, "w")
    out.write("Error\tCount\tFraction\n")
    for tag in ordered_tags:
        count = _error_count.get(tag, 0)
        fraction = _error_count.get(tag, 0)/sum(_error_count.values())
        out.write("{}\t{}\t{}\n".format(tag, count, fraction))
    out.close()

    file_name = summary_path + ".cell_barcodes.csv"
    sorted_barcodes = sorted(_cell_count, key=lambda x: _cell_count[x],
            reverse=True)
    cell_cumsum = cumsum([_cell_count[b] for b in sorted_barcodes]) / \
            sum(_cell_count.values())
    out = open(file_name, "w")
    out.write("Cell_Barcode\tCount\tCumulative_Sum\n")
    for i, barcode in enumerate(sorted_barcodes):
        out.write("{}\t{}\t{}\n".format(barcode, _cell_count[barcode],
            cell_cumsum[i]))
    out.close()

def main():
    args = parse_args()

    global _tag_bc
    global _tag_umi
    global _tag_error
    _tag_bc    = args.tag_bc
    _tag_umi   = args.tag_umi
    _tag_error = args.tag_error

    bam_write_mode = "w" if args.output == "-" else "wb"
    in_filename1, in_filename2 = args.input

    # Processing ####
    _start = time.perf_counter()
    in_seqs1 = (seq for head, seq, qual in FastqGeneralIterator(gzopen(in_filename1, "rt")))
    in_reads2 = (record for record in FastqGeneralIterator(gzopen(in_filename2, "rt")))

    if args.subset: # DEBUGGING PURPOSES
        in_seqs1 = islice(in_seqs1, args.subset)

    if args.pipeline == "dropseq":
        bam_header = {'HD':{'VN': '1.6', 'SO':'unknown'}}
        out_bam = pysam.AlignmentFile(args.output, bam_write_mode, header=bam_header)

    elif args.pipeline == "scpipe":
        out_fastq = gzopen(args.output, "wt")

    logging.info("Extracting tags.")
    pool = Pool(args.cores)
    for (i, tags) in enumerate(pool.imap(get_tags, in_seqs1), 1):
        if (i) % 1e6 == 0:
            logging.info("{} reads processed.".format(str(i)))
        if args.pipeline == "dropseq":
            title, seq, qual = next(in_reads2)
            sam_record = pysam.AlignedSegment()
            sam_record.query_name = title.split()[0]
            sam_record.query_sequence = seq
            sam_record.query_qualities = pysam.qualitystring_to_array(qual)
            sam_record.template_length = len(seq)
            sam_record.flag = 4
            sam_record.set_tags(tags.items())
            out_bam.write(sam_record)
        elif args.pipeline == "scpipe" and tags.get(_tag_bc):
            title, seq, qual = next(in_reads2)
            out_fastq.write("@{}_{}#{}\n{}\n+\n{}\n".format(tags[_tag_bc], tags[_tag_umi], title.split()[0], seq, qual))

        if args.summary_prefix: # summary statistics
            compute_summary(tags)
    logging.info("{} reads processed.".format(str(i)))
    pool.close()
    logging.info("All reads analyzed.")

    if args.pipeline == "dropseq":
        out_bam.close()
    elif args.pipeline == "scpipe":
        out_fastq.close()

    if args.summary_prefix:
        logging.info("Writing summary files.")
        write_summary(args.summary_prefix)

    _end = time.perf_counter()
    logging.info("Done. Elapsed time: {} minutes.".format(round((_end-_start)/60, 2)))

def parse_args():
    description = "A tool to extract cellular and molecular identifiers " +\
        "from single cell RNA sequencing experiments"

    parser = ArgumentParser(description=description)

    parser.add_argument("-i", "--input", nargs=2,
            help="R1 and R2 from a paired end sequencing experiment")
    parser.add_argument("-o", "--output", required=True,
            help="Tagged unmapped BAM file (use '-' to output to stdout")

    parser.add_argument("--pipeline", default="dropseq", choices=["dropseq", "scpipe"],
        help="Set output type depending on pipeline tool chosen")

    parser.add_argument("-s", "--summary-prefix",
        help="Prefix for summary files (including absolute or relative paths)")

    parser.add_argument("-c","--cores", type=int, default=1,
        help="Number of processing units (CPUs) to use (default=1)")

    parser.add_argument("--tag-bc", type=str, default="XC",
        help="Tag for single cell barcode (default=XC)")

    parser.add_argument("--tag-umi", type=str, default="XM",
        help="Tag for Unique Molecular Identifier (default=XM)")

    parser.add_argument("--tag-error", type=str, default="XE",
        help="Tag for errors (default=XE)")

    parser.add_argument("--subset", type=int,
        help="Select a lower number of reads to analyze [debugging]")

    parser.add_argument("-v", '--version', action='version', version='%(prog)s 0.9.0')

    args = parser.parse_args()

    global _barcodes
    barcodes_file = Path(sys.argv[0]).resolve().parent.joinpath("barcodes.txt")
    _barcodes = barcodes_file.open().read().split()

    # check parameters
    if not (args.tag_bc != args.tag_umi != args.tag_error != args.tag_bc):
        logging.error("Tags provided with '--tag' must be all different.")
        sys.exit(1)

    tag_pattern = "^[A-Za-z][A-Za-z0-90]$"
    if not (re_match(tag_pattern, args.tag_bc) and re_match(tag_pattern, args.tag_umi) and re_match(tag_pattern, args.tag_error)):
        logging.error("Tags provided with '--tag' must be two-character strings matching /[A-Za-z][A-Za-z0-9]/")
        sys.exit(1)

    for x in args.input:
        ext = "".join(Path(x).suffixes)
        if ext != ".fastq.gz":
            logging.error("Invalid input file extensions '{}': '.fastq.gz' is required".format(ext))
            sys.exit(1)

    return(args)

if __name__ == "__main__":
    main()
