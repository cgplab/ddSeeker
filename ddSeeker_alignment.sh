#!/usr/bin/env bash
# MIT License
#
# Copyright 2017 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

set -eu
# Fail if any of the commands in a pipeline fails
set -o pipefail

tmpdir=
outdir=`pwd`
genomedir=
reference=
pipeline=0
ncores=1
echo_prefix=
dropseq_root=$(dirname $0)
star_executable=STAR
progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform tagging, filtering and alignment

-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: directory containing this script.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-g <genomedir>      : Directory of STAR genome directory.  Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-s <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
-k <ddSeeker_path>  : Full path of ddSeeker.  Default: ddSeeker is found via PATH environment variable.
-n <num-cores>      : Number of CPUs unit for parallel computation.  Default: 1.
-t <tmpdir>         : Where to write temporary files.  Default: a new subdirectory in /tmp
-e                  : Echo commands instead of executing them.  Cannot use with -p.
-p                  : Reduce file I/O by pipeline commands together.  Requires more memory and processing power.
-h                  : Print usage and exit.
EOF
}

function error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]; then
        error_exit "$name has not been specified.  $flag flag is required"
    fi
}


while getopts ":d:k:g:n:o:r:s:t:eph" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    o ) outdir=$OPTARG;;
    g ) genomedir=$OPTARG;;
    r ) reference=$OPTARG;;
    k ) ddseeker_executable=$OPTARG;;
    s ) star_executable=$OPTARG;;
    n ) ncores=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    e ) echo_prefix="echo";;
    p ) pipeline=1;;
    h ) usage
            exit 1;;
    \?) usage
            exit 1;;
    * ) usage
            exit 1;;
  esac
done
shift $(($OPTIND - 1))

if [[ "$pipeline" == 1 && -n "$echo_prefix" ]]; then
    error_exit "-p and -e cannot be used together"
fi

check_set "$dropseq_root" "Drop-seq root" "-d"
check_set "$genomedir" "Genome directory" "-g"
check_set "$reference" "Reference fasta"  "-r"

if (( $# != 1 )); then
    error_exit "Incorrect number of arguments"
fi

if [[ -z "$tmpdir" ]]; then
    tmpdir=`mktemp -d`
    echo "Using temporary directory $tmpdir"
fi

if [[ "$star_executable" != "STAR" ]]; then
    if [[ ! ( -x $star_executable && -f $star_executable ) ]]; then
        error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
    fi
elif which STAR > /dev/null; then
    echo > /dev/null
else
    error_exit "STAR executable must be on the path"
fi

refflat=${reference%.fa*}.refFlat
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar

unmapped_bam=$1
unaligned_tagged_bam=${tmpdir}/unaligned_tagged.bam
unaligned_tagged_filtered_bam=${tmpdir}/unaligned_tagged_filtered.bam
unaligned_tagged_filtered_fastq=${tmpdir}/unaligned_tagged_filtered.fastq
aligned_sam=${tmpdir}/star.Aligned.out.sam
aligned_sorted_bam=${tmpdir}/aligned_sorted.bam
files_to_delete="${unaligned_tagged_filtered_bam} ${aligned_sorted_bam} ${aligned_sam}"

# Stage 1: pre-alignment tag and trim
tag_reads="${ddseeker_executable} ${unmapped_bam} --summary ${outdir}/unaligned_tagged_summary --ncores ${ncores}"
filter_bam="${dropseq_root}/FilterBAM TAG_REJECT=XE"

# Stage 2: alignment
sam_to_fastq="java -Xmx500m -jar ${picard_jar} SamToFastq
    I=${unaligned_tagged_filtered_bam}"
star_align="$star_executable --genomeDir ${genomedir} --runThreadN $ncores
    --outFileNamePrefix ${tmpdir}/star."

# Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
sort_aligned="java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50
    -XX:GCHeapFreeLimit=10 -Xmx4000m -jar ${picard_jar} SortSam
    I=${aligned_sam} O=${aligned_sorted_bam} SORT_ORDER=queryname
    TMP_DIR=${tmpdir}"

# Stage 4: merge and tag aligned reads
merge_bam="java -Xmx4000m -jar ${picard_jar} MergeBamAlignment
    REFERENCE_SEQUENCE=${reference} UNMAPPED_BAM=${unaligned_tagged_filtered_bam}
    ALIGNED_BAM=${aligned_sorted_bam} INCLUDE_SECONDARY_ALIGNMENTS=false
    PAIRED_RUN=false"
tag_with_gene_exon="${dropseq_root}/TagReadWithGeneExon
    O=${outdir}/star_gene_exon_tagged.bam ANNOTATIONS_FILE=${refflat} TAG=GE"

if (( $pipeline == 1 )); then
    # Stage 1
    $tag_reads --output-bam "-" | \
        $filter_bam I=/dev/stdin O=${unaligned_tagged_filtered_bam}

    # Stage 2
    $sam_to_fastq FASTQ=/dev/stdout | \
        $star_align --readFilesIn /dev/stdin

    # Stage 3
    $sort_aligned

    # Stage 4
    $merge_bam O=/dev/stdout COMPRESSION_LEVEL=0 | \
        $tag_with_gene_exon I=/dev/stdin
else
    # Stage 1
    $echo_prefix $tag_reads --output-bam $unaligned_tagged_bam
    $echo_prefix $filter_bam I=$unaligned_tagged_bam O=$unaligned_tagged_filtered_bam
    files_to_delete="$files_to_delete $unaligned_tagged_bam"

    # Stage 2
    $echo_prefix $sam_to_fastq FASTQ=$unaligned_tagged_filtered_fastq
    $echo_prefix $star_align --readFilesIn $unaligned_tagged_filtered_fastq
    files_to_delete="$files_to_delete $unaligned_tagged_filtered_fastq"

    # Stage 3
    $echo_prefix $sort_aligned

    # Stage 4
    $echo_prefix $merge_bam O=$tmpdir/merged.bam
    $echo_prefix $tag_with_gene_exon I=$tmpdir/merged.bam
    files_to_delete="$files_to_delete $tmpdir/merged.bam"

fi

$echo_prefix rm $files_to_delete
