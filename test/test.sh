#!/bin/bash
set -eu

../ddSeeker.py unmapped.bam \
    -o tagged.bam \
    -s summary \
    -n 10
