#!/usr/bin/env python3

# Date: 2019-03-22
# Author: Dario Romagnoli
# Script to split cell line bam files into single cell bamfiles

import simplesam
import time
from argparse import ArgumentParser
from os import mkdir
from os.path import exists
from os.path import join as join_path
from sys import argv, exit, stderr

def info(message):
    current_time = time.localtime()
    print('[{time}] {message}'.format(
        time=time.strftime('%Y-%m-%d %H:%M:%S', current_time),
        message=message), file=stderr)

def pretty(d):
    for key, value in d.items():
        print("* {key}: {value}".format(key=str(key), value=str(value)), file=stderr)

def main():
    parser = ArgumentParser(description = "Utility to split a tagged BAM/SAM " + \
            "file into separate SAM files for different barcodes.")
    parser.add_argument("input", help = "BAM/SAM tagged file")
    parser.add_argument("outdir", help = "Directory where to store the splitted SAM files")
    parser.add_argument("-t", "--tag-name", default="XC", help = "TAG name (default=XC)")
    parser.add_argument("-b", "--barcodes-file", help = "Selection of cell barcodes (default: all barcodes in BAM/SAM file)")
    args = parser.parse_args()

    mkdir(args.outdir)
    # user provided a subset of barcodes (create dictionary and set a flag)
    if args.barcodes_file:
        barcodes = dict.fromkeys(open(args.barcodes_file).read().split(), 0)
        subset_barcodes = True
    else:
        barcodes = dict()
        subset_barcodes = False

    in_samfile = open(args.input, "r")
    in_sam = simplesam.Reader(in_samfile)

    print("Selected TAG name: {name}".format(name=args.tag_name), file=stderr)
    info("Analyzing file {file}".format(file=args.input))

    selected_reads = 0
    for tot_reads, read in enumerate(in_sam):
        if (tot_reads + 1) % 1e5 == 0:
            info("Reads (selected/total): {s}/{t}".format(s=str(selected_reads+1), t=str(tot_reads+1)))

        cell_bc = read.tags.get(args.tag_name)
        if (cell_bc and subset_barcodes and (cell_bc in barcodes)) or (cell_bc and not subset_barcodes):
            selected_reads += 1
            out_filename = join_path(args.outdir, cell_bc + ".sam")
            if not barcodes.get(cell_bc): # NB: both get zero and get null equal false
                barcodes[cell_bc] = barcodes.get(cell_bc, 0) + 1
                with open(out_filename, "w") as out_samfile:
                    out_sam = simplesam.Writer(out_samfile, in_sam.header)
                    out_sam.write(read)
            else:
                barcodes[cell_bc] += 1
                with open(out_filename, "a") as out_samfile:
                    print(str(read), end="", file=out_samfile)
    in_sam.close()

if __name__ == "__main__":
    main()
