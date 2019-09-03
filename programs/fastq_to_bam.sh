#!/bin/bash

set -e

pg=`basename ${BASH_SOURCE[0]}`

# get the size of memory and cpu cores
function mem_size
{
	if [[ $(command -v free) ]]; then
		free | grep -i '^Mem' | awk '{print $2}'
	else
		exit 7
	fi
}

function num_cores
{
	if [[ $(command -v nproc) ]]; then
		echo `nproc`
	else
		exit 8
	fi
}

# test whether a bam file is paired end
function paired_bam
{
    if [[ $( samtools view -f 0x1 $1) ]]; then
            echo "ok"
    else
            echo "";
    fi
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

# check whether a bam is duplicate-marked
function dup_marked
{
	cnt=$(samtools view $1 | head -10000 | grep 'PG:Z:MarkDuplicates' | wc -l)
	if [[ $cnt -gt 0 ]]; then
		echo "D";
	else
		echo "";
	fi
}

# remove duplicate reads from a bam file
function remove_dups()
{
	bam=$1
	markedBam=${bam/%.bam/.dup_marked.bam}
	noDupBam=${bam/%.bam/.no_dup.bam}
	metric=${bam/%.bam/.dup_metrics.txt}
	# if dup marked already, just filter
	if [[ $(dup_marked $bam) ]]; then
		warn "$bam is already dup-marked"
		samtools view -F 0x400 -b -o $noDupBam $bam
		echo "$noDupBam"
		return
	fi
	# otherwise mark dups
	if [[ ! $(bam_sorted $bam) ]]; then
		samtools sort -m $memPerCore -@ $extraCores -o tmp.$$.bam $bam
		mv tmp.$$.bam $bam
		msg "$bam is coordinate-sorted now"
	fi
	$picard MarkDuplicates \
	I=$bam ASO=coordinate M=$metric O=$markedBam \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 TMP_DIR=. \
	TAGGING_POLICY=All
	mv $markedBam $bam; # change original bam to marked one
	msg "$bam is dup-marked now"
	samtools view -F 0x400 -b -o $noDupBam $bam
	echo "$noDupBam"
}

function usage
{
	cat <<EOF
Usage: $0 [options] <fastq-file1> [<fastq-file2>]

This program accepts 1 or 2 fastq files and does the following steps:

1. trim adaptors and low-quality bases from input fastq files.
2. align trimmed fastq files to a reference genome.
3. remove duplicates (see option --remove-dups) and convert 
the alignment bams to bed format.

If only 1 fastq file is provided, then the program is run in
single-end mode. If two fastq files, then run in pair-end mode.

Options (default values are in []):

-h|--help:  show this usage message.

--sample-id: sample id, used for naming all output files.
    [string matching /^[a-zA-Z0-9_]+/ in input filename ]

--trim-galore: the path to program 'trim_galore'. [trim_galore]

--q-cutoff: minimal base quality value for base trimming. [10]

--min-len: minimal length for trimmed reads to keep. [20]

--trim-out: output folder for trimmed sequence files. [trimmed]

--tg-params: other parameters directly fed to trim-galore. []

--bowtie: path to the program 'bowtie2'. [bowtie2]

--n-threads: number of threads to run bowtie. Consider change
default value when you run multiple resource-demanding programs
at the same time [$alnCores].

--ref-index: the filepath to bowtie index files, just as for the
option '-x' for bowtie2. []

--bowtie-params: other parameters passed to bowtie, such as 
'-k 4 -X 700 --local', must be quoted. []

--remove-dups: if provided, PCR/picard duplicates are marked in
the generated bam file and will be skipped when converting bam
file to bed format. Default is keeping duplicate reads.

--bam-filter: further parameters to filter bam files during the
conversion to bed format. These parameters are fed to 'samtools',
so they need be in quotes like '-q 30'. Check 'samtools view' for
available filtering options. []

--resume: from which step the program starts to run. Available values
are 'trim', 'align', 'to-bed'. This option allows the program to pick
up from where they were stopped. When specifying this option, all the
other options should be exactly the same as last run, otherwise errors
would occur [trim]

Examples:
# run the program on a pair of read files
$0 --trim-galore mydir/trim_galore --ref-index bowtie2_index/hg38 --q-cutoff 20 --trim-out . fastq/test_1.fastq.gz fastq/test_2.fastq.gz

# run the program on a pair of read files with output filenames
# starting with 'mysample'
$0 --sample-id mysample --trim-galore mydir/trim_galore --ref-index bowtie2_index/hg38 --q-cutoff 20 --trim-out . fastq/test_1.fastq.gz fastq/test_2.fastq.gz

# the save as above, but give parameters to bowtie
$0 --trim-galore mydir/trim_galore --bowtie2-params '--local -X2000'  --ref-index bowtie2_index/hg38 --q-cutoff 20 --trim-out . fastq/test_1.fastq.gz fastq/test_2.fastq.gz

# resume the above run from step 'align'
$0 --trim-galore mydir/trim_galore --ref-index bowtie2_index/hg38 --q-cutoff 20 --trim-out . --resume align fastq/test_1.fastq.gz fastq/test_2.fastq.gz

# run the program on a pair of read files and remove duplicates when
# converting to bed files
$0 --trim-galore mydir/trim_galore --ref-index bowtie2_index/hg38 \
--q-cutoff 20 --trim-out . --remove-dups fastq/test_1.fastq.gz fastq/test_2.fastq.gz

# ***** ENCODE ATAC-seq parameters ***** #
$0 --trim-galore mydir/trim_galore --ref-index bowtie2_index/hg38 \
--q-cutoff 20 --trim-out . --remove-dups -bowtie-params \
'-X2000 -k 4 --local' --bam-filter '-F 1804 -q 30 -f 2' --tg-params \
'-e 0.2' --min-len 5 fastq/test_1.fastq.gz fastq/test_2.fastq.gz

EOF

}

function msg
{
	echo "[$pg] $*" >&2
}

function warn
{
	msg "[WARNING: $pg] $*"
}

function clean_up
{
	if [[ $tmpFiles ]]; then
		rm -rf ${tmpFiles[@]}
	fi
}

tmpFiles=();

## get system resources
mem=$(( `mem_size`*1024 )); # RAM size in bytes
usableMem=`echo "scale=0; $mem*0.8" | bc`
nCores=`num_cores`;
alnCores=$nCores; # cores for bowtie

if [[ $# -lt 1 ]]; then
	usage;
	exit 1;
fi

msg "Input arguments: $*"

# global variables and their default values
paired=""; # reads are in paired mode?
sampleId=''
bamDir=bams
beginStep='trim'
## trim parameters
tg='trim_galore'
trimOut='trimmed'
baseQualCut=10; # base quality cutoff for trim
minReadLen=20
tgParams="";
## align parameters
bowtie='bowtie2'
bowtieParams="";
refBowtieIndex=''
## convert to bed file
removeDups=''
bamFilter=''

# read arguments
while [[ $# -gt 0 ]]
do
	k=$1; # key
	shift; # remove the key from the stack
	case $k in 
		-h|--help)
			usage;
			exit 1;
			;;
		--trim-galore)
			tg=$1;
			shift;
			;;
		--sample-id)
			sampleId=$1;
			shift;
			;;
		--q-cutoff)
			baseQualCut=$1;
			shift;
			;;
		--min-len)
			minReadLen=$1;
			shift;
			;;
		--trim-out)
			trimOut=$1;
			shift;
			;;
		--tg-params)
			tgParams=$1;
			shift;
			;;
		--bowtie)
			bowtie=$1;
			shift;
			;;
		--n-threads)
			alnCores=$1;
			shift;
			;;
		--ref-index)
			refBowtieIndex=$1;
			shift;
			;;
		--bowtie-params)
			bowtieParams=$1;
			shift;
			;;
		--remove-dups)
			removeDups="T";
			;;
		--bam-filter)
			bamFilter=$1;
			shift;
			;;
		--resume)
			beginStep=$1;
			shift;
			;;
		*) # unknown parameters
			posArgs+=("$k") # save arguments to array
			;;
	esac

done

depends=($tg cutadapt $bowtie samtools bedtools)

for e in ${depends[@]}
do
	if [[ ! $(command -v $e) ]]; then
		msg "Can't find command $e"
		exit 2;
	fi
done

# find the program picard if duplicates need be marked
if [[ "$removeDups" == "T" ]]; then
	if [[ $(command -v picard) ]]; then
		picard='picard'
	elif [[ $(command -v picard-tools) ]]; then
		picard="picard-tools"
	else
		warn "Can't find 'picard' or 'picard-tools'"
		exit 2;
	fi
fi

if [[ ! "$refBowtieIndex" ]]; then
	msg "Bowtie index of reference genome is needed"
	exit 2;
fi

if [[ ${#posArgs[@]} -lt 1 ]]; then
	msg "No input fastq files"
	exit 3
fi

if [[ ${#posArgs[@]} -gt 2 ]]; then
	msg "At most 2 fastq files can be accepted"
	msg "But the provided parameters are '${posArgs[@]}'"
	exit 4;
fi

fq1=${posArgs[0]}
fq2=${posArgs[1]}

if [[ "$fq2" ]]; then
	paired="T"
	msg "This program is run in paired-ended mode"
else
	paired="";
	msg "This program is run in single-ended mode"
fi

if [[ ! "$sampleId" ]]; then
	# get default sample Id
	f=$(basename $fq1)
	[[ $f =~ ^([a-zA-Z0-9_]+) ]];
	sampleId=${BASH_REMATCH[1]}
	msg "Set sample id to '$sampleId'"
fi

stepNum=0;

case $beginStep in
	trim)
		stepNum=1;
		;;
	align)
		stepNum=2;
		;;
	to-bed)
		stepNum=3
		;;
	*)
		msg "Unknown resume step '$beginStep'"
		exit 6;
		;;
esac

msg "Job started at `date`"
extraCores=$(( alnCores - 1 )); # cores for samtools
memPerCore=$(echo "scale=0; $usableMem/$nCores" | bc);
msg "System resources: #cores $nCores; memory, $mem"
msg "Used: #cores: $alnCores; memory per core, $memPerCore"

mkdir -p $bamDir

msg "Step 1: trimming fastq files at `date`"

#paired="--paired"
mode=""
if [[ "$paired" == "T" ]]; then
	mode="--paired"
	fqs="$fq1 $fq2"
else
	fqs="$fq1"
fi

mkdir -p $trimOut
msg "*** trim-galore parameters: $tgParams"
if [[ $stepNum -le 1 ]]; then
	$tg -q $baseQualCut --gzip --length $minReadLen --trim-n \
	$mode -o $trimOut --basename $sampleId $tgParams $fqs
else
	msg "Skipping 'trim' step in resume mode"
fi

# run the trimming now
#cmd="$tg ${posArgs[@]}"
#$tg -q $baseQualCut
#
#msg "running '$cmd'"
#$tg ${posArgs[@]}

msg "Step 2: align reads to reference genome $refBowtieIndex at `date`"

if [[ "$paired" ]]; then
	fq1=$trimOut/${sampleId}_1_val_1.fq.gz;
	fq2=$trimOut/${sampleId}_2_val_2.fq.gz;
	if [[ ! -f $fq1 ]] || [[ ! -f $fq2 ]]; then
		warn "Trimmed fastq files '$fq1' and '$fq2' can't be found"
		exit 5;
	fi
	fqs="-1 $fq1 -2 $fq2"
else
	msg "Single-end alignment has not implemented yet"
	exit 5;
fi

commOpts="$bowtieParams -p $alnCores"

#msg "**** Bowtie common parameters: [$localAln] $commOpts"
#msg "*** bam filter '$bamFilter'"

refBam=$bamDir/$sampleId.bam

if [[ $stepNum -le 2 ]]; then
	$bowtie $commOpts -x $refBowtieIndex $fqs \
	| samtools view -b -o $refBam
	msg "Reference bam file is generated: $refBam"
else
	msg "Skipping 'align' step in resume mode"
fi

msg "Step 3: convert bam files to bed files at `date`"

if [[ $stepNum -le 3 ]]; then
	if [[ $removeDups ]]; then
		msg "Mark and skip duplicates from '$refBam'"
		refBam=$(remove_dups $refBam)
		tmpFiles+=($refBam)
	fi
	for bam in "$refBam"
	do
		bedFile=$sampleId.aligned.bed.gz
	
		if [[ "$paired" == "T" ]]; then
			# sort bam file according to names
			tmpBam=tmp.$$.bam
			samtools sort -n -o $tmpBam -@ $extraCores -m $memPerCore $bam
			samtools view -u $bamFilter $tmpBam | bamToBed -bedpe | \
			gawk 'BEGIN{OFS="\t"}$1==$4{print $1,$2,$6,$7,$8,$9 "/" $10}' \
			| sort -k1,1 -k2,2n -k3,3n | gzip -c >$bedFile
			rm $tmpBam
		else
			samtools view -u $bamFilter $bam | bamToBed \
			| sort -k1,1 -k2,2n -k3,3n | gzip -c >$bedFile
		fi
		
		msg "The bam-to-bed file '$bedFile' is generated"
	done
else
	msg "Skipping 'to-bed' step in resume mode"
fi

# remove temporary files
clean_up;

msg "Job done at `date`"

exit 0;

