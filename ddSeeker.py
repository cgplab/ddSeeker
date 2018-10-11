#!/usr/bin/env python3

import time
import argparse
from re import match
import sys
from os.path import abspath, dirname
from os.path import join as pjoin
from numpy import cumsum
from multiprocessing import Pool
from itertools import islice
import itertools
import pysam
from Bio import pairwise2

local_alignment = pairwise2.align.localxs
global_alignment = pairwise2.align.globalxs
_linkers = ["TAGCCATCGCATTGC", "TACCTCTGAGCTGAA"]

def info(*args, **kwargs):
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), file=sys.stderr, *args, **kwargs)

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
        score = global_alignment(bc, block, -1, -1,
                                 score_only=True,
                                 one_alignment_only=True)
        if len(block) >= 5 and score >= 5:
            return(bc)
    else:
        return(None)

def get_tags(sequence):
    """Extract barcodes from R1 and return a tuple of SAM-format TAGs.
    XB = barcode
    XU = UMI
    XE = error

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
        alignment = local_alignment(sequence, linker, -2, -1,
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
        return([(tag_error, "LX")]) # no linker aligned
    elif not starts[0]:
        return([(tag_error, "L1")]) # linker 1 not aligned
    elif not starts[1]:
        return([(tag_error, "L2")]) # linker 2 not aligned

    if starts[1]-starts[0] == 21+k[0]:
        bc2 = sequence[starts[1]-6: starts[1]]
    elif starts[1]-starts[0] == 20+k[0]: # 1 deletion in bc2
        bc2 = sequence[starts[1]-5: starts[1]]
    elif starts[1]-starts[0] == 22+k[0]: # 1 insertion in bc2
        bc2 = sequence[starts[1]-7: starts[1]]
    else:
        return([(tag_error, "I")])

    if starts[0] < 5:
        return([(tag_error, "D")])
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
        return([(tag_error, "J")])

    gac = sequence[starts[1]+32+k[1]: starts[1]+35+k[1]]
    try:
        dist_gac = hamming_dist(gac, "GAC")
    except ValueError:
        dist_gac = float("inf")
    if dist_gac > 1:
        return([(tag_error, "K")])

    bc3 = sequence[starts[1]+15+k[1]: starts[1]+21+k[1]]

    barcode = []
    for block in (bc1, bc2, bc3):
        fixed = fix_block(block)
        if fixed:
            barcode.append(fixed)
        else:
            return([(tag_error, "B")])

    umi = sequence[starts[1]+24+k[1]: starts[1]+32+k[1]]

    return([(tag_bc, "".join(barcode)), (tag_umi, umi)])

def add_tags(paired_reads):
    read1, read2 = paired_reads
    tags = get_tags(read1.seq)
    for tag in tags:
        read2[tag[0]] = tag[1]
    read2.flag -= (1 + 8 + 128) # remove these flags: read paired, mate unmapped, second in pair
    return(read2)

def main(args):
    args = parse_args(args)
    in_filename = args.input_bam
    out_filename = args.output_bam
    summary = args.summary
    cores = args.ncores
    barcodes_file = args.barcodes_file

    tag_pattern = "^[A-Za-z][A-Za-z0-90]$"
    if not (args.tag_bc != args.tag_umi != args.tag_error != args.tag_bc):
        print("Error: tags must be unique.", file=sys.stderr)
        sys.exit(1)
    if not (match(tag_pattern, args.tag_bc) and match(tag_pattern, args.tag_umi) and match(tag_pattern, args.tag_error)):
        print("Error: tags must be two-character strings matching /[A-Za-z][A-Za-z0-9]/", file=sys.stderr)
        sys.exit(1)

    global tag_bc
    global tag_umi
    global tag_error
    tag_bc    = args.tag_bc
    tag_umi   = args.tag_umi
    tag_error = args.tag_error

    global _barcodes
    try:
        _barcodes = [_.rstrip().split()[0] for _ in open(barcodes_file).readlines()[:96]]
    except FileNotFoundError:
        print("Error: '{}' file not found.".format(barcodes_file) + \
             "Specify file path with -b flag or run 'ddSeeker_barcodes.py' to create one.", file=sys.stderr)
        sys.exit(1)

    #  info("Start analysis:", in_filename, "->", "stdout" if out_filename == "-" else out_filename)
    in_bam1 = pysam.AlignmentFile(in_filename, "rb", check_sq=False)
    in_bam2 = pysam.AlignmentFile(in_filename, "rb", check_sq=False)
    out_bam = pysam.AlignmentFile(out_filename, "w" if out_filename == "-" else "wb", template = in_bam1)

    # iterable of read1 sequences
    in_bam_iter1 = islice(in_bam1.fetch(until_eof=True), None, None, 2)
    in_bam_iter2 = islice(in_bam2.fetch(until_eof=True), 1, None, 2)
    reads = (_.query_sequence for _ in in_bam_iter1)

    p = Pool(cores)
    #  info("Get identifiers from R1")
    cell_count = {}
    error_count = {}
    _start = time.perf_counter()
    for tags in p.imap(get_tags, reads):
        read2 = next(in_bam_iter2)
        read2.set_tags(tags)
        read2.flag -= (1 + 8 + 128) # remove flags: read paired, mate unmapped, second in pair
        out_bam.write(read2)

        # summary statistics
        if summary and tags[0][0] == tag_error:
            error_count[tags[0][1]] = error_count.get(tags[0][1], 0) + 1
        elif summary and tags[0][0] == tag_bc:
            cell_count[tags[0][1]] = cell_count.get(tags[0][1], 0) + 1
            error_count["PASS"] = error_count.get("PASS", 0) + 1

    _check = time.perf_counter()
    print("Finished writing bam file in {} seconds.".format(round(_check-_start, 3)))
    in_bam1.close()
    in_bam2.close()
    out_bam.close()
    p.close()

    if summary:
        info("Write summary files")
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
    _end = time.perf_counter()
    print("Finished writing summary file in {} seconds.".format(round(_end-_check, 3)))

def parse_args(args):
    """Parse argv"""
    description = "A tool to extract cellular and molecular identifiers " +\
        "from single cell RNA sequencing experiments"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("input_bam", help="Merged paired-end uBAM file")

    parser.add_argument("-o", "--output-bam", default="-",
        help="Tagged uBAM file (default=<stout>")

    parser.add_argument("-b", "--barcodes-file",
        default=pjoin(dirname(abspath(sys.argv[0])), "barcodes.txt"),
        help="Barcode blocks file (default=<ddSeeker_path>/barcodes.txt")

    parser.add_argument("-s", "--summary",
        help="Summary files name prefix (including absolute or relative path)")

    parser.add_argument("-n","--ncores", type=int, default=1,
        help="Number of processing units (CPUs) to use (default=1)")

    parser.add_argument("--tag-bc", type=str, default="XB",
        help="Tag for single cell barcode (default=XB)")

    parser.add_argument("--tag-umi", type=str, default="XU",
        help="Tag for Universal Molecular Identifier (default=XU)")

    parser.add_argument("--tag-error", type=str, default="XE",
        help="Tag for errors (default=XE)")

    parser.add_argument("-v", '--version', action='version', version='%(prog)s 1.2')

    args = parser.parse_args()
    return(args)

if __name__ == "__main__":
    main(sys.argv[1:])
