#!/usr/bin/env python3

import logging
import sys
import time
from argparse import ArgumentParser, ArgumentTypeError
import pysam
from pathlib import Path
from functools import partial
from multiprocessing import Pool
from itertools import islice
from Bio import pairwise2, SeqIO
from gzip import open as gzopen

logging.basicConfig(level=logging.INFO, datefmt='%H:%M:%S',
                    format="[%(asctime)s] %(levelname)s - %(message)s")
_match_table = [1, 0, -1, -1]
_linkers = ["TAGCCATCGCATTGC", "TACCTCTGAGCTGAA"]
_local_aligner = partial(pairwise2.align.localms, one_alignment_only=True)

class countDict(dict):
    """
    >>> x = countDict([("a",1), ("b",1)])
    >>> y = countDict([("a",1), ("b",2), ('c',10)])
    >>> x+y
    {'a': 2, 'c': 10, 'b': 3}

    >>> x = countDict([("a",1), ("b",2)])
    >>> print(x)
    {'a':1, 'b':2}
    >>> x.update(["a", "b", "c"])
    >>> print(x)
    {'a':2, 'c':1, 'b':3}

    """
    def update(self, iterable):
        for key in iterable:
            self[key] = self.get(key, 0) + 1

def extract_barcodes(read):
    starts = []
    for linker in _linkers:  # align the two linkers
        alignment = _local_aligner(read, linker, *_match_table)[0]
        check_length = alignment[4] - alignment[3] == len(linker)  # end - start
        check_match = alignment[2] >= len(linker) - 1  # alignment score (allow 1 mismatch)
        if (check_length and check_match):
            starts.append(alignment[3])
        else:
            starts.append("")

    if all(starts) and (starts[1]-starts[0] == 21) and starts[0] >= 6:
        bc_blocks = [str(read[starts[0] - 6: starts[0]]),
                     str(read[starts[1] - 6: starts[1]]),
                     str(read[starts[1] + 15: starts[1] + 21])]
        if any([len(_) != 6 for _ in bc_blocks]):
            bc_blocks = []
    else:
        bc_blocks = []

    return(bc_blocks)

def main(argv=None):
    if not argv:
        argv = sys.argv[1:]

    # Parameters ####
    args = parse_args(argv)
    in_files = args.input
    out_file = args.output
    cores = args.cores

    bc_count = countDict()

    # Processing ####
    _start = time.perf_counter()
    for in_file in in_files:
        ext = Path(in_file).suffixes
        if ext[-1] == ".bam" or ext[-1] == ".sam":
            file_type = "bam"
        elif ext[-2] == ".fastq" and ext[-1] == ".gz":
            file_type = "fastq"
        else:
            raise IOError("Unrecognized input file extension '{}'".format(ext[-1]))

        logging.info("Extracting barcodes from '{}'".format(in_file))
        if file_type == "bam":
            in_bam = pysam.AlignmentFile(in_file, check_sq=False)
            in_bam_iter = islice(in_bam.fetch(until_eof=True), None, None, 2)
            reads = (_.query_sequence.upper() for _ in in_bam_iter)
        else:
            reads = (str(_.seq.upper()) for _ in SeqIO.parse(gzopen(in_file, "rt"), format="fastq"))

        if args.subset: # DEBUGGING PURPOSES
            reads = islice(reads, args.subset)

        pool = Pool(cores)
        barcodes = pool.map(extract_barcodes, reads)
        for bc_blocks in barcodes:
            bc_count.update(bc_blocks)
        pool.close()

        if file_type == "bam":
            in_bam.close()

    sorted_bc = sorted(bc_count, key=lambda x: bc_count[x], reverse=True)
    with out_file.open("w") as out:
        for k in sorted_bc:
            out.write("{}\t{}\n".format(k, bc_count[k]))
    _end = time.perf_counter()
    logging.info("Done. Elapsed time: {} minutes.".format(round((_end-_start)/60, 2)))

def parse_args(args):
    description = "A tool to generate a set of barcode blocks " +\
            "to be used by ddSeeker.py"

    parser = ArgumentParser(description=description)

    parser.add_argument("-i", "--input", required=True, nargs="+",
            help = "Merged paired-end uBAM file(s) or fastq files")

    barcodes_file = Path(sys.argv[0]).resolve().parent.joinpath("barcodes.txt")
    parser.add_argument("-o", "--output", type=Path, default=barcodes_file,
            help = "File reporting barcodes from most represented to less represented" +\
                   "(default=/path-to-ddSeeker/barcode.txt)")

    parser.add_argument("-c","--cores", type=int, default=1,
            help="Number of processing units (CPUs) to use (default=1)")

    parser.add_argument("--subset", type=int,
        help="Select a lower number of reads to analyze [for debugging purposes]")

    parser.add_argument("-v", '--version', action='version', version='%(prog)s 0.9.0')
    args = parser.parse_args()

    return(args)

if __name__ == "__main__":
    main()
