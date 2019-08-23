# This file contains common NGS functions and can be sourced into
# other bash files.

## testing bam files
function is_bam_paired
{
        if [[ $( samtools view -f 0x1 $1 | head -2) ]]; then
                echo "ok"
        else
                echo "";
        fi
}

function is_bam_sorted
{
        line=$(samtools view -H $1 | grep '^@HD' | gawk 'BEGIN{FS="\t"}$3=="SO:coordinate"')
        if [[ $line ]]; then
                echo "S";
        else
                echo "";
        fi
}

## summarizing bams
# count the number of read pairs as well reads in a bam file
function count_read_pairs
{
        bam=$1
        # check whether all reads in paired mode, even only one end is
        # mapped
        unpairedCnt=$(samtools view -c -F 0x1 $bam)
        if [[ $unpairedCnt -gt 0 ]]; then
                msg "There are $unpairedCnt unpaired reads in $bam. can't proceed"
                exit 4;
        fi
        unmappedCnt=$(samtools view -c -f 0x4 $bam)
        if [[ $unmappedCnt -gt 0 ]]; then
                msg "There are $unmappedCnt unmapped reads in $bam. can't proceed"
                exit 4;
        fi
        # get all reads (each name counts only once)
        nReads=$(samtools view -c -F 0x100 $bam)
        single=$(samtools view -c -F 0x100 -f 0x8 $bam)
        nPairs=$(echo "scale=0; ($nReads-$single)/2+$single" | bc -l) # each pair counted only once
        echo "$nPairs $nReads"
}

