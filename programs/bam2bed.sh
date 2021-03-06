#!/bin/bash

function usage
{
        cat <<EOF
Usage: $0 <bam-file> [<out-file>]

This programs read a bam file and convert it to a bed file
using bedtools. When the <bam-file> contains any paired-end
reads (based on flag 0x1), the whole bam file is processed
in paired-end mode.

The out-file name is optional. Default is print to standard
output.

E.g.: $0 test.bam | gzip -c >test.bed.gz

EOF
}


function msg
{
        echo "$*" >&2
}

# test whether a bam file is paired end
function paired_bam
{
        if [[ $( samtools view -f 0x1 $1 | head -2) ]]; then
                echo "ok"
        else
                echo "";
        fi
}

if [[ $# -lt 1 ]]; then
        usage
        exit 1;
fi

bam=$1
out=$2

if [[ $(paired_bam $bam) ]]; then
        msg "Converting $bam to bed in paired-end mode"
        # need sort the bam file first
        tmpBam=tmp.$$.bam
        samtools sort -n -o $tmpBam $bam
        bam=$tmpBam
        if [[ $out ]]; then
                samtools view -F 0x4 -u $bam | bamToBed -bedpe -i - | gawk 'BEGIN{OFS="\t"}$1==$4{print $1,$2,$6,$7,$8,$9 "/" $10}' >$out
        else
                samtools view -F 0x4 -u $bam | bamToBed -bedpe -i - | gawk 'BEGIN{OFS="\t"}$1==$4{print $1,$2,$6,$7,$8,$9 "/" $10}'
        fi
        rm $tmpBam
else
        msg "Converting $bam to bed in single-end mode"
        if [[ $out ]]; then
                samtools view -F 0x4 -u $bam | bamToBed -i - >$out
        else
                samtools view -F 0x4 -u $bam | bamToBed -i -
        fi
fi

msg "Job is done"

exit 0;
