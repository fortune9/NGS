#!/bin/bash

nCPUs=`nproc`
outDir="dup_marked"

function usage()
{

    cat <<EOF
Usage: $0 [options] <bam-file1> [<bam-file2> [...]]

This program marks duplicated reads in input bam files.
For each bam file, it will generate two files: duplicate-
marked bam file and metrics file.

Options (default values are in []):

--out-dir <path>: the directory to store output files. [$outDir]

--n-cpus <int>: the number of cpus to use if the program
        'parallel' is available [$nCPUs]

Example use:
$0 --out-dir test.bam

EOF
}

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

if [[ ! $(command -v picard) ]]; then
    echo "Command 'picard' is required"
    exit 3
fi

bamFiles=()

while [[ $# -gt 0 ]];
do
    k=$1; shift;
    case $k in
        --out-dir)
            outDir=$1;
            shift;
            ;;
        --n-cpus)
            nCPUs=$1;
            shift;
            ;;
        *)
            bamFiles+=("$k")
            ;;
    esac
done

if [[ ${#bamFiles[@]} -lt 1 ]]; then
    echo "Bam files are required"
    exit 2;
fi

if [[ ${#bamFiles[@]} -lt $nCPUs ]]; then
    nCPUs=${#bamFiles[@]}
fi

mkdir -p $outDir

if [[ $(command -v parallel) ]]; then
    echo "Processing files with gnu parallel at `date`"
    parallel -j $nCPUs picard MarkDuplicates \
        I={} O=$outDir/{/.}.dupMarked.bam M=$outDir/{/.}.dupMetrics.txt \
        TAGGING_POLICY=All TMP_DIR=. ::: ${bamFiles[@]}
else
    echo "Processing files with for loop at `date`"
    for f in ${bamFiles[@]}
    do
        outBam=`basename ${f/%.bam/.dupMarked.bam}`
        outMetrics=`basename ${f/%.bam/.dupMetrics.txt}`

        picard MarkDuplicates \
            I=$f O=$outBam M=$outMetrics \
            TAGGING_POLICY=All TMP_DIR=.
    done
fi

echo "Job done at `date`"

exit 0

