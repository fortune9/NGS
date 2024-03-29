#!/bin/bash

set -e

#tg=/home/ec2-user/tools/extra/tmp-work/TrimGalore-0.6.1/trim_galore
#bowtie_index=/home/ec2-user/tools/extra/tmp-work/bowtie2_index/hg38.no_alt/hg38

# this program accepts a bam file and run the following steps
# 1. subsample 5m reads and estimate library size
# 2. count the number of reads mapped to chrM and chrY
# 3. call peaks on the subsampled bam file

function usage
{
	cat <<EOF
Usage: $0 [options] <bamFile>  

Examples:
# sample 10% of input reads
$0 -o test_out -f 0.1 --paired test.bam
# sample 100 reads in a paired-end mode, i.e., 50 read pairs
E.g.: $0 -o test_out -f 100 --paired test.bam
# the same as above, but read from AWS S3
E.g.: $0 -o test_out -f 100 --paired s3://bucket1/to/test.bam

This program reads a bam file, remove duplicates, filter properly
aligned read pairs (paired-end mode only), and sample a specified
number of reads from it. Optionally, it can also remove reads from
chromosomes chrM and chrY.

Input bam file can be local or at AWS S3.

options:

-h/--help:  show this help message.
--paired:   reads in the bam are in paired-end. Default is single-end.
-f/--frac:  the fraction of reads sampled from input bam. If greater
            than 1, the parameter is assumed as the absolute number of reads
            to sample; in paired-end mode, if one wants to sample
            n read pairs, he need specify 2*n. 1 means no sampling.  [1].
-d/--keep-dup: if provided, reads from PCR/optical duplicates will be
            retained
-b/--bad-aligned: if provided, unproperly aligned paired-end reads
            will be retained.
--no-MY:    if provided, reads on chromosome chrM and chrY will be
            removed.
--mapq:     a number, reads with MAPQ smaller than this will be
            removed before sampling. Default is 0.
-t/--keep-temp: if provided, the intermediate files will be kept after
            run.
-o/--out:   the basename for the output. Default is the input bam
            filename with suffix '.bam' replaced with ".sub"
--seed:     the random number seed for sampling. Default is current
			PID.

EOF

}

function msg
{
	echo -e "$*" >&2;
}

# get the file if they are located remotely
# parameters: input file, the variable name to store new file
function get_file() {
    f=$1
    if [[ $f =~ ^s3://. ]]; then
        newF=$(basename $f)
        echo "Downloading $f"
        aws s3 cp --quiet $f $newF && \
            echo "Downloading $f is successful"
        s3dir=$(dirname $f)
        eval "$2=$newF"
        tmpFiles+=($newF)
    fi
}

# upload file to s3
# parameters: the name of a file to upload
function upload_file() {
    aws s3 cp --quiet $1 $s3dir/$1 && \
        echo "Uploading $1 to $s3dir/$1 is successful"
}

# check whether a bam is coordinate sorted
function bam_sorted
{
	line=$(samtools view -H $1 | grep '^@HD' | gawk 'BEGIN{FS="\t"}$3=="SO:coordinate"')
	if [[ $line ]]; then
		echo "S";
	else
		echo "";
	fi
}

# count the number of read pairs as well reads in a bam file
function count_read_pairs
{
	bam=$1
    nCpus=${2:-0}
	# check whether all reads in paired mode, even only one end is
	# mapped
	unpairedCnt=$(samtools view -c -F 0x1 -@ $nCpus $bam)
	if [[ $unpairedCnt -gt 0 ]]; then
		msg "There are $unpairedCnt unpaired reads in $bam. can't proceed"
		exit 4;
	fi
	unmappedCnt=$(samtools view -c -f 0x4 -@ $nCpus $bam)
	if [[ $unmappedCnt -gt 0 ]]; then
		msg "There are $unmappedCnt unmapped reads in $bam. can't proceed"
		exit 4;
	fi
	# get all reads (each name counts only once)
	nReads=$(samtools view -c -F 0x100 -@ $nCpus $bam)
	single=$(samtools view -c -F 0x100 -f 0x8 -@ $nCpus $bam)
	nPairs=$(echo "scale=0; ($nReads-$single)/2+$single" | bc -l) # each pair counted only once
	echo "$nPairs $nReads"
}

# check whether a command exists
function check_exe
{
	if [[ $(command -v $1) ]]; then
		echo "ok"
	else
		echo "";
	fi
}

# convert read/pair number to fraction
function num_to_frac
{
	bam=$1;
	num=$2;
	total=$( samtools idxstats $bam | gawk 'BEGIN{s=0}{s+=$3}END{print s}' )
	echo "scale=$(( ${#total} + 1 )); $num/$total" | bc
}

function cmp_num()
{
	if [[ $(echo "$1 $3 $2" | bc) -gt 0 ]]; then
		echo 1
	else
		echo ""
	fi
}

function delete_files() {
    for f in "$@"
    do
        if [[ -e "$f" ]]; then
            echo Deleting "$f"
            rm "$f"
        else
            echo "$f doesn't exist"
        fi
    done
}

if [[ $# -lt 1 ]]; then
	usage;
	exit 1;
fi

if [[ $(check_exe picard) ]]; then
	picard='picard'
elif [[ $(check_exe picard-tools) ]]; then
	picard="picard-tools"
else
	msg "can't find picard or picard-tools in the system"
	exit 2;
fi

depends=(samtools)
for i in ${depends[@]}
do
	ok=$(check_exe $i);
	if [[ ! $ok ]]; then
		msg "command $i can't be found"
		exit 2;
	fi
done

## process command arguments
posArgs=();
frac=1;
keepDup="";
badAlign="";
keepTmp="";
outBase="";
noMY="";
mapqCut=0;
paired="";

while [[ $# -gt 0 ]]
do
	k=$1;
	shift; # remove this value from stack
	#echo $k "->" $1 
	case $k in
		-h|--help)
			usage;
			exit 1;
			;;
		--paired)
			paired=T;
			;;
		-f|--frac)
			frac="$1"
			shift; # remove the value
			;;
		-d|--keep-dup)
			keepDup="T"
			;;
		-b|--bad-aligned)
			badAlign=T
			;;
		--no-MY)
			noMY=T
			;;
		--mapq)
			mapqCut=$1;
			shift;
			;;
		-t|--keep-temp)
			keepTmp=T
			;;
		-o|--out)
			outBase="$1";
			shift;
			;;
		*) # unknown options
			posArgs+=("$k") # save in array for later use
			;;
	esac
done

set -- "${posArgs[@]}"
bam=$1;
s3dir=""
tmpFiles=();
get_file $bam bam
## globals
ncores=`nproc`
extraCpus=$(( ncores - 1))
memPerCore=2G
seed=$$;

# set default outbase
if [[ ! $outBase ]]; then
	outBase=$(basename ${bam/%.bam/.sub}).$$
fi

mkdir -p $(dirname $outBase)

# if input bam is not sorted, sort it first
if [[ ! $(bam_sorted $bam) ]]; then
	samtools sort -m $memPerCore -@ $extraCpus -o tmp.$$.bam $bam
	mv tmp.$$.bam $bam
	msg "$bam is coordinate-sorted now"
fi

if [[ ! $keepDup ]]; then
	msg "Marking and removing PCR duplicates"
	markedBam=${outBase}.dup_marked.bam
	metric=${outBase}.dup_metrics.txt
	$picard MarkDuplicates \
	I=$bam ASO=coordinate M=$metric O=$markedBam \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 TMP_DIR=. \
	TAGGING_POLICY=All
	#aws s3 cp --quiet $markedBam $s3dirOut/$markedBam
	#aws s3 cp --quiet $metric $s3dirOut/$metric
	bam=$markedBam
	tmpFiles+=($markedBam);
fi
# count the number of reads before filtering
readCntOrig=($(count_read_pairs $bam $extraCpus)); # store in an array
echo "[STAT] In original file: there are ${readCntOrig[1]} reads in ${readCntOrig[0]} pairs"
# index the bam file
samtools index -@ $extraCpus $bam
tmpFiles+=($bam.bai)
filterParams="";

if [[ $mapqCut -gt 0 ]]; then
	msg "Reads with MAPQ < $mapqCut to be removed"
	filterParams+=" -q $mapqCut"
fi

if [[ (! $badAlign) && $paired ]]; then
	msg "Unproperly aligned read pairs to be removed"
	filterParams+=" -f 0x2"
fi

if [[ $noMY ]]; then
	msg "Reads on chrM/chrY to be removed"
	samtools idxstats $bam | grep -vP '^(chrM|chrY)' \
	| gawk 'BEGIN{OFS="\t"}{print $1,0,$2,$1}' >tmp1.$$.bed
	filterParams+=" -L tmp1.$$.bed"
	tmpFiles+=(tmp1.$$.bed);
fi

if [[ ! $keepDup ]]; then # remove duplicates
	filterParams+=" -F 0x400"
fi

# filter reads now
if [[ $filterParams ]]; then
	msg "Filter bam with options '$filterParams' "
	o=tmp.$$.filtered.bam
	samtools view $filterParams -@ $extraCpus -b -o $o $bam
	bam=$o
	samtools index $o
	tmpFiles+=($o $o.bai)
fi

# count reads after filtering
readCntFiltered=($(count_read_pairs $bam $extraCpus)); # store in an array
echo "[STAT] In filtered file: there are ${readCntFiltered[1]} reads in ${readCntFiltered[0]} pairs"
# now sub-sample the reads

if [[ $(cmp_num $frac 1 ">") -gt 0 ]]; then
	msg "Converting reads number to fraction of the total"
	frac=$(num_to_frac $bam $frac)
fi

if [[ $(cmp_num $frac 1 ">=") -gt 0 ]]; then
	msg "Fraction $frac is equal/greater than 1, so no sampling will be done"
	exit 2;
fi

frac=${frac/#*\.} # remove anything until decimal point
if [[ ! "$frac" =~ ^[0-9]+$ ]]; then
	msg "Fraction digits '$frac' contains non-numbers"
	exit 3;
fi

subfrac=$seed.$frac
subBam=$outBase.bam
echo "[STAT] Sampling 0.$frac fraction reads from filtered bam"
samtools view -@ $extraCpus -b -o $subBam -s $subfrac $bam
msg "Final bam file is generated: $subBam"

# count reads after filtering
readCntSampled=($(count_read_pairs $subBam $extraCpus)); # store in an array
echo "[STAT] In sampled file: there are ${readCntSampled[1]} reads in ${readCntSampled[0]} pairs"

ratio1=$(echo "scale=4; ${readCntFiltered[1]}/${readCntOrig[1]}" | bc -l)
ratio2=$(echo "scale=4; ${readCntFiltered[0]}/${readCntOrig[0]}" | bc -l) # pairs ratio
echo -e "The filtered reads/pairs ratio (against the original): $ratio1\t$ratio2"
ratio1=$(echo "scale=4; ${readCntSampled[1]}/${readCntOrig[1]}" | bc -l)
ratio2=$(echo "scale=4; ${readCntSampled[0]}/${readCntOrig[0]}" | bc -l) # pairs ratio
echo -e "The sampled reads/pairs ratio (against the original): $ratio1\t$ratio2"

if [[ ! $keepTmp ]] && [[ ${#tmpFiles[@]} -gt 0 ]]; then # clean up
	msg "Removing temporary files"
    delete_files "${tmpFiles[@]}"
fi

if [[ "$s3dir" ]]; then
    upload_file $subBam
    # remove the middle files
    if [[ ! $keepTmp ]]; then
        delete_files $subBam
    fi
fi

msg "Job is done at `date`"

exit 0;

