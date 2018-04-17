#!/usr/bin/python3

import time
import argparse
from sys import argv
from os.path import abspath, dirname
from os.path import join as pjoin
from multiprocessing import Pool
from itertools import islice, chain
import pysam
from Bio import pairwise2
from Bio.Seq import Seq

local_alignment = pairwise2.align.localms
_linkers = ["TAGCCATCGCATTGC", "TACCTCTGAGCTGAA"]
_match_table = [1, 0, -1, -1]

def info(*args):
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), *args)

def extract_barcodes(read):
    """
    """
    read = Seq(read)
    starts = []
    for linker in _linkers:  # align the two linkers
        alignment = local_alignment(read, linker, *_match_table,
                one_alignment_only=True)[0]
        check_length = alignment[4] - alignment[3] == len(linker)  # end - start
        check_match = alignment[2] >= len(linker) - 1  # alignment score (allow 1 mismatch)
        if (check_length and check_match):
            starts.append(alignment[3])
        else:
            starts.append("")

    if all(starts) and (starts[1]-starts[0] == 21) and starts[0] >= 6:
        bc1 = str(read[starts[0] - 6: starts[0]])
        bc2 = str(read[starts[1] - 6: starts[1]])
        bc3 = str(read[starts[1] + 15: starts[1] + 21])
        if len(bc1) == 6 and len(bc2) == 6 and len(bc3) == 6:
            return([bc1, bc2, bc3])
        else:
            return([])
    else:
        return([])

def main(args):
    args = parse_args(args)
    in_files = args.input_files
    out_file = args.output_file
    cores = args.ncores

    bc_count = {}
    for in_file in in_files:
        info("Extract barcodes from", in_file)
        in_bam = pysam.AlignmentFile(in_file, "rb", check_sq=False)
        in_bam_iter = islice(in_bam.fetch(until_eof=True), None, None, 2)
        reads = (_.query_sequence for _ in in_bam_iter)
        # parallel processing of reads
        pool = Pool(cores)
        barcodes = pool.map(extract_barcodes, reads)
        pool.close()
        in_bam.close()
        info("Count barcodes")
        for bc in chain.from_iterable(barcodes):
            bc_count[bc] = bc_count.get(bc, 0) + 1
        del(barcodes)

    sorted_bc = sorted(bc_count, key=lambda x: bc_count[x], reverse=True)
    info("Write output")
    out = open(out_file, "w")
    i = 0
    while i < len(sorted_bc):
        k = sorted_bc[i]
        out.write("{}\t{}\n".format(k, bc_count[k]))
        i += 1
    out.close()
    info("Done")

def parse_args(args):
    """Parse argv"""
    description = "A tool to generate a set of barcode blocks " +\
            "to be used by ddSeeker.py"
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("input_files", nargs="+",
            help = "Merged paired-end uBAM files")
    barcodes_file = pjoin(dirname(abspath(argv[0])), "barcodes.txt")
    parser.add_argument("-o", "--output", dest="output_file", default=barcodes_file,
            help = "File reporting barcodes from most represented to less represented" +\
                   "(default=/path-to-ddSeeker/barcode.txt)")
    parser.add_argument("-n","--ncores", type=int, default=1,
            help="Number of processing units (CPUs) to use (default=1)")
    parser.add_argument("-v", '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()
    return(args)

if __name__ == "__main__":
    main(argv[1:])
