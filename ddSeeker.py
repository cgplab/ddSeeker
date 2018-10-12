#!/usr/bin/env python3

import logging
import sys
import time
import argparse
import pysam
from re import match as re_match
from os.path import abspath, dirname
from os.path import join as pjoin
from numpy import cumsum
from multiprocessing import Pool
from itertools import islice
from Bio import pairwise2, SeqIO
from gzip import open as gzopen
from itertools import islice

logging.basicConfig(level=logging.INFO, datefmt='%H:%M:%S',
                    format="[%(asctime)s] %(levelname)s - %(message)s")
_local_aligner = pairwise2.align.localxs
_global_aligner = pairwise2.align.globalxs
_linkers = ["TAGCCATCGCATTGC", "TACCTCTGAGCTGAA"]

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
        score = _global_aligner(bc, block, -1, -1,
                                 score_only=True,
                                 one_alignment_only=True)
        if len(block) >= 5 and score >= 5:
            return(bc)
    else:
        return(None)

def get_tags(sequence):
    """Extract barcodes from R1 and return a tuple of SAM-format TAGs.
    Default tag values:
    - XB = barcode
    - XU = UMI
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
        alignment = _local_aligner(sequence, linker, -2, -1,
                one_alignment_only=True)[0]
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
        return([(_tag_error, "LX")]) # no linker aligned
    elif not starts[0]:
        return([(_tag_error, "L1")]) # linker 1 not aligned
    elif not starts[1]:
        return([(_tag_error, "L2")]) # linker 2 not aligned

    if starts[1]-starts[0] == 21+k[0]:
        bc2 = sequence[starts[1]-6: starts[1]]
    elif starts[1]-starts[0] == 20+k[0]: # 1 deletion in bc2
        bc2 = sequence[starts[1]-5: starts[1]]
    elif starts[1]-starts[0] == 22+k[0]: # 1 insertion in bc2
        bc2 = sequence[starts[1]-7: starts[1]]
    else:
        return([(_tag_error, "I")])

    if starts[0] < 5:
        return([(_tag_error, "D")])
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
            return([("XE", "J")])

    try:
        dist_acg = hamming_dist(acg, "ACG")
    except ValueError:
        dist_acg = float("inf")
    if dist_acg > 1:
        return([(_tag_error, "J")])

    gac = sequence[starts[1]+32+k[1]: starts[1]+35+k[1]]
    try:
        dist_gac = hamming_dist(gac, "GAC")
    except ValueError:
        dist_gac = float("inf")
    if dist_gac > 1:
        return([(_tag_error, "K")])

    bc3 = sequence[starts[1]+15+k[1]: starts[1]+21+k[1]]

    barcode = []
    for block in (bc1, bc2, bc3):
        fixed = fix_block(block)
        if fixed:
            barcode.append(fixed)
        else:
            return([(_tag_error, "B")])

    umi = sequence[starts[1]+24+k[1]: starts[1]+32+k[1]]

    return([(_tag_bc, "".join(barcode)), (_tag_umi, umi)])

def compute_summary(tags):
    # summary statistics
    if tags[0][0] == _tag_error:
        error_count[tags[0][1]] = error_count.get(tags[0][1], 0) + 1
    elif tags[0][0] == _tag_bc:
        cell_count[tags[0][1]] = cell_count.get(tags[0][1], 0) + 1
        error_count["PASS"] = error_count.get("PASS", 0) + 1

def write_summary(summary):
    file_name = summary + ".errors"
    ordered_tags = ["LX", "L1", "L2", "I", "D", "J", "K", "B", "PASS"]
    out = open(file_name, "w")
    out.write("Error\tCount\tFraction\n")
    for tag in ordered_tags:
        count = error_count.get(tag, 0)
        fraction = error_count.get(tag, 0)/sum(error_count.values())
        out.write("{}\t{}\t{}\n".format(tag, count, fraction))
    out.close()

    file_name = summary + ".cell_barcodes"
    sorted_barcodes = sorted(cell_count, key=lambda x: cell_count[x],
            reverse=True)
    cell_cumsum = cumsum([cell_count[b] for b in sorted_barcodes]) / \
            sum(cell_count.values())
    out = open(file_name, "w")
    out.write("Cell_Barcode\tCount\tCumulative_Sum\n")
    for i, barcode in enumerate(sorted_barcodes):
        out.write("{}\t{}\t{}\n".format(barcode, cell_count[barcode],
            cell_cumsum[i]))
    out.close()

def main(argv=None):
    if not argv:
        argv = sys.argv[1:]

    # Parameters #####
    args = parse_args(argv)

    if len(args.input) == 1:
       in_filename1, in_filename2 = args.input[0], None
    elif len(args.input) == 2:
        in_filename1, in_filename2 = args.input
    else:
        print("no more than two inputs.")
        sys.exit(1)

    out_filename = args.output
    summary_pref = args.summary_prefix
    cores = args.cores
    barcodes_file = args.barcodes_file

    global _tag_bc
    global _tag_umi
    global _tag_error
    _tag_bc    = args.tag_bc
    _tag_umi   = args.tag_umi
    _tag_error = args.tag_error

    global _barcodes
    try:
        _barcodes = [_.rstrip().split()[0] for _ in open(barcodes_file).readlines()[:96]]
    except FileNotFoundError:
        print("'{}' file not found.".format(barcodes_file) + \
              "  Specify file path with -b flag or run 'ddSeeker_barcodes.py' to create one.", file=sys.stderr)
        sys.exit(1)

    if summary_pref:
        global cell_count
        global error_count
        cell_count = {}
        error_count = {}

    # Processing ####
    _start = time.perf_counter()
    pool = Pool(cores)
    bam_write_mode = "w" if out_filename == "-" else "wb"
    if not in_filename2: # unmapped bam file
        in_bam1 = pysam.AlignmentFile(in_filename1, "rb", check_sq=False)
        in_bam_iter1 = islice(in_bam1.fetch(until_eof=True), None, None, 2)
        reads = (_.query_sequence for _ in in_bam_iter1)

        in_bam2 = pysam.AlignmentFile(in_filename1, "rb", check_sq=False)
        in_bam_iter2 = islice(in_bam2.fetch(until_eof=True), 1, None, 2)

        out_bam = pysam.AlignmentFile(out_filename, bam_write_mode, template=in_bam1)
    else: # fastq files
        reads = (str(_.seq) for _ in SeqIO.parse(gzopen(in_filename1, "rt"), format="fastq"))
        in_fastq2 = (_ for _ in SeqIO.parse(gzopen(in_filename2, "rt"), format="fastq"))

        header = {'HD':{'VN': '1.6', 'SO':'unknown'}}
        out_bam = pysam.AlignmentFile(out_filename, bam_write_mode, header=header)

    if args.subset: # DEBUGGING PURPOSES
        reads = islice(reads, args.subset)

    logging.info("Extracting tags")
    for tags in pool.imap(get_tags, reads):
        if not in_filename2:
            read2 = next(in_bam_iter2)
            sam_record = read2
        else:
            read2 = next(in_fastq2)
            sam_record = pysam.AlignedSegment()
            sam_record.query_name = read2.name
            sam_record.query_sequence = str(read2.seq)
            sam_record.query_qualities = pysam.qualitystring_to_array(read2.format("fastq").split()[-1])
        #  sam_record.flag -= (1 + 8 + 128) # remove flags: read paired, mate unmapped, second in pair
        sam_record.flag = 4
        sam_record.template_length=len(read2.seq)
        sam_record.set_tags(tags)
        out_bam.write(sam_record)

        if summary_pref: # summary statistics
            compute_summary(tags)
    out_bam.close()
    logging.info("All reads analyzed")

    if not in_filename2:
        in_bam1.close()
        in_bam2.close()
    out_bam.close()
    pool.close()

    if summary_pref:
        logging.info("Writing summary files")
        write_summary(summary_pref)

    _end = time.perf_counter()
    logging.info("Done. Elapsed time: {} minutes.".format(round((_end-_start)/60, 2)))

def parse_args(args):
    description = "A tool to extract cellular and molecular identifiers " +\
        "from single cell RNA sequencing experiments"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", "--input", required=True, nargs="*",
            help="Either one merged paired-end unmapped BAM file " +\
                 "or two paired-end FASTQ files")
    parser.add_argument("-o", "--output", required=True,
            help="Tagged unmapped BAM file (use '-' to output to stdout")

    parser.add_argument("-b", "--barcodes-file",
        default=pjoin(dirname(abspath(sys.argv[0])), "barcodes.txt"),
        help="Barcode blocks file (default=<ddSeeker_path>/barcodes.txt")

    parser.add_argument("-s", "--summary-prefix",
        help="Prefix for summary files (including absolute or relative paths)")

    parser.add_argument("-c","--cores", type=int, default=1,
        help="Number of processing units (CPUs) to use (default=1)")

    parser.add_argument("--tag-bc", type=str, default="XB",
        help="Tag for single cell barcode (default=XB)")

    parser.add_argument("--tag-umi", type=str, default="XU",
        help="Tag for Universal Molecular Identifier (default=XU)")

    parser.add_argument("--tag-error", type=str, default="XE",
        help="Tag for errors (default=XE)")

    parser.add_argument("--subset", type=int,
        help="Select a lower number of reads to analyze [for debugging purposes]")

    parser.add_argument("-v", '--version', action='version', version='%(prog)s 0.90.0')

    args = parser.parse_args()

    # check parameters
    if not (1 <= len(args.input) <= 2):
        raise argparse.ArgumentTypeError("'-i/--input' expects exactly 1 or 2 files")

    if not (args.tag_bc != args.tag_umi != args.tag_error != args.tag_bc):
        raise argparse.ArgumentTypeError("tags provided with '--tag' must be unique.")

    tag_pattern = "^[A-Za-z][A-Za-z0-90]$"
    if not (re_match(tag_pattern, args.tag_bc) and re_match(tag_pattern, args.tag_umi) and re_match(tag_pattern, args.tag_error)):
        raise argparse.ArgumentTypeError("tags provided with '--tag' must be two-character strings matching /[A-Za-z][A-Za-z0-9]/")
    return(args)

if __name__ == "__main__":
    main()
