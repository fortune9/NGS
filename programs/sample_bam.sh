#!/bin/bash

set -e

s3dirIn=s3://zymo-filesystem/home/zzhang/work-partition/ATAC-seq/Analysis/in1408
s3dirOut=s3://zymo-filesystem/home/zzhang/work-partition/ATAC-seq/Analysis/in1408
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

E.g.: $0 -o test_out -f 0.1 --paired test.bam # 10% reads
E.g.: $0 -o test_out -f 100 --paired test.bam # 100 reads, or 50 pairs

This program reads a bam file, remove duplicates, filter properly
aligned read pairs (paired-end mode only), and sample a specified
number of reads from it. Optionally, it can also remove chrM and chrY
reads.

For paired-end bam, if option --frac is fed with a number greater than
1, then this number is regarded as the number of reads, not read
pairs, so a number 1000 means 500 read pairs.

options:

-h/--help:  show this help message.
--paired:   reads in the bam are in paired-end. Default is single-end.
-f/--frac:  the fraction of reads sampled from input bam. If greater
            than 1, the parameter is assumed as the absolute number of reads
            required. 1 means no sampling [1].
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
	echo "$*" >&2;
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

if [[ $# -lt 1 ]]; then
	usage;
	exit 1;
fi

depends=(samtools picard)
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
## globals
ncores=`nproc`
extraCpus=$(( ncores - 1))
memPerCore=2G
tmpFiles=();
seed=$$;

# set default outbase
if [[ ! $outBase ]]; then
	outBase=$(basename ${bam/%.bam/.sub})
fi

if [[ ! $keepDup ]]; then
	msg "Marking and removing PCR duplicates"
	markedBam=${outBase}.dup_marked.bam
	metric=${outBase}.dup_metrics.txt
	if [[ ! $(bam_sorted $bam) ]]; then
		samtools sort -m $memPerCore -@ $extraCpus -o tmp.$$.bam $bam
		mv tmp.$$.bam $bam
		msg "$bam is coordinate-sorted now"
	fi
	picard MarkDuplicates \
	I=$bam ASO=coordinate M=$metric O=$markedBam \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 TMP_DIR=. \
	TAGGING_POLICY=All
	#aws s3 cp --quiet $markedBam $s3dirOut/$markedBam
	#aws s3 cp --quiet $metric $s3dirOut/$metric
	bam=$markedBam
	tmpFiles+=($markedBam $markedBam.bai);
fi

# index the bam file
samtools index $bam
filterParams="";

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
	samtools view -@ $extraCpus -b $filterParams -q $mapqCut -o $o $bam
	bam=$o
	samtools index $o
	tmpFiles+=($o $o.bai)
fi

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
msg "Sampling 0.$frac fraction reads from filtered bam"
samtools view -@ $extraCpus -b -o $subBam -s $subfrac $bam

if [[ ! $keepTmp ]]; then # clean up
	msg "Removing temporary files"
	rm "${tmpFiles[@]}"
fi

msg "Job is done at `date`"

exit 0;

	metric=$alnDir/$id.dup_metrics.$sizeTag.$i.txt
	samtools view -@ $extraCpus -b -o $subBam -s $subfrac $bam
	samtools index $subBam
	# estimate library size
	picard EstimateLibraryComplexity I=$subBam O=$metric TMP_DIR=.
	# count reads
	chrMYCnt=$alnDir/${id}.reads_cnt.$sizeTag.$i.tsv
	samtools idxstats $subBam >tmp1.$$.txt
	echo -e "sample\tmappedCnt.all\tmappedCnt.chrM\tmappedCnt.chrY" >$chrMYCnt
	chrMcnt=$( less tmp1.$$.txt | gawk 'BEGIN{s=0}$1~/^chrM[^\d]?/{s+=$3}END{print s}' )
	chrYcnt=$( less tmp1.$$.txt | gawk 'BEGIN{s=0}$1~/^chrY[^\d]?/{s+=$3}END{print s}' )
	allcnt=$( less tmp1.$$.txt | gawk 'BEGIN{s=0}{s+=$3}END{print s}' )
	echo -e "$id\t$allcnt\t$chrMcnt\t$chrYcnt" >>$chrMYCnt

size=5000000
sizeTag=5m
id=`basename ${bam%.bam}`
if [[ $outBase ]]; then
	id=$outBase
fi
s3dir=""
species=hs

alnDir=sub_aln.$$
macsDir=sub_macs.$$
mkdir -p $alnDir $macsDir

echo Job started at `date`

# preprocess bams
if [[ "$dupMarked" == "T" ]]; then
	echo The input has been dup-marked, so no more mark.
else
	echo "Step 0: Mark duplicates on $bam"
	markedBam=$alnDir/${id}.dup_marked.bam
	metric=$alnDir/${id}.dup_metrics.txt
	samtools sort -m 2G -@ $extraCpus -o tmp.$$.bam $bam
	mv tmp.$$.bam $bam
	picard MarkDuplicates \
	I=$bam ASO=coordinate M=$metric O=$markedBam \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 TMP_DIR=. \
	TAGGING_POLICY=All
	#aws s3 cp --quiet $markedBam $s3dirOut/$markedBam
	#aws s3 cp --quiet $metric $s3dirOut/$metric
	bam=$markedBam
fi

# start sub-sampling and analyzing

# get the fraction to sample
echo Subsampling $frac fraction from $bam
frac=${frac/#*\.}

for i in `seq $numReps`
do
	echo "Step 1 (rep $i): sub-sample $size reads from $bam"
	# use $i as seed
	subfrac=$i.$frac
	subBam=$alnDir/$id.$sizeTag.$i.bam
	metric=$alnDir/$id.dup_metrics.$sizeTag.$i.txt
	samtools view -@ $extraCpus -b -o $subBam -s $subfrac $bam
	samtools index $subBam
	# estimate library size
	picard EstimateLibraryComplexity I=$subBam O=$metric TMP_DIR=.
	# count reads
	chrMYCnt=$alnDir/${id}.reads_cnt.$sizeTag.$i.tsv
	samtools idxstats $subBam >tmp1.$$.txt
	echo -e "sample\tmappedCnt.all\tmappedCnt.chrM\tmappedCnt.chrY" >$chrMYCnt
	chrMcnt=$( less tmp1.$$.txt | gawk 'BEGIN{s=0}$1~/^chrM[^\d]?/{s+=$3}END{print s}' )
	chrYcnt=$( less tmp1.$$.txt | gawk 'BEGIN{s=0}$1~/^chrY[^\d]?/{s+=$3}END{print s}' )
	allcnt=$( less tmp1.$$.txt | gawk 'BEGIN{s=0}{s+=$3}END{print s}' )
	echo -e "$id\t$allcnt\t$chrMcnt\t$chrYcnt" >>$chrMYCnt
	## copy files
	#aws s3 cp --quiet $subBam $s3dirOut/$subBam
	#aws s3 cp --quiet $metric $s3dirOut/$metric
	#aws s3 cp --quiet $chrMYCnt $s3dirOut/$chrMYCnt
	
	echo "Step 2 (rep $i): call peaks on $subBam"
	
	# filter the bams
	if [[ $paired ]]; then
		filter="-f 0x2"
		format=BAMPE
	else
		filter="";
		format=AUTO
	fi
	o=tmp.$$.filtered.bam
	samtools idxstats $subBam | grep -vP '^(chrM|chrY)' \
	| gawk 'BEGIN{OFS="\t"}{print $1,0,$2,$1}' >tmp1.$$.bed
	samtools view -@ $extraCpus -b $filter -F 0x400 -q 30 -L tmp1.$$.bed -o $o $subBam
	macs2 callpeak -t $o -f $format -g $species --outdir $macsDir -B \
	-n "$id.$sizeTag.$i" --trackline -q 0.05 --call-summits \
	--keep-dup all
done

#aws s3 sync --quiet macs $s3dirOut/macs

# delete all generated files
#rm -rf $alnDir $macsDir
rm tmp1.$$.txt tmp1.$$.bed tmp.$$.filtered.bam

#rm tmp1.$$.txt $subBam $metric $chrMYCnt

echo "Results are put in folders $alnDir and $macsDir"

echo Job is done at `date`;

exit 0;

