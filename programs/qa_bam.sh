#!/bin/bash

set -e

#pg=`basename ${SOURCE_BASH[0]}`
tmpFiles=();

function rand_str
{
	len=${1:-8}; # string length
	echo $(head -10 /dev/urandom | tr -dc '0-9a-zA-Z' | fold -w $len | head -1)
}

## let's load common functions
function download_file
{
	out=tmp.$$.$(rand_str)
	wget -O $out --quiet "$1" || (msg "downloding '$1' failed "; fail 3);
	echo $out
}

bashFuncUrl='https://raw.githubusercontent.com/fortune9/programs/master/bash_functions.sh'
ngsFuncUrl='https://raw.githubusercontent.com/fortune9/NGS/master/ngs_functions.sh'
bashFuncFile=$(download_file $bashFuncUrl)
ngsFuncFile=$(download_file $ngsFuncUrl)
source $bashFuncFile
source $ngsFuncFile
tmpFiles+=($bashFuncFile $ngsFuncFile)

function usage
{
	cat <<EOF
Usage: $0 [options] <bam-file>

This program does a quality analysis of the input bam file, including
the following steps:

1. markdup: generate a report using picard's MarkDuplicates.
2. flagstat: output a flagstat report using samtools.
3. libcomplex: compute library complexity (NRF, PBC1, PBC2).

Options (default values are in []):

--skip <string>: what steps to skip. The steps can be specified using
both numbers and names, such as '1,2' and 'markdup,flagstat' are
equivalent. []

--out <string>: output filename. The default is writing to screen [].

--bam-filter <string>: the parameters used by samtools to filter the
input bam before running step 2, such as '-f 0x2 -F 0x700', must
be quoted, and check 'samtools view' for other filter options. []

** libcomplex options **

--chrs <string>: in addition to the metrics for the bam file overall,
the total and unique reads for these chromosomes will also be
separately output. The format should be like 'chr1,chr2'. []

Example uses:

$0 test.bam

## skip step 3
$0 --skip 3 test.bam

## filter input bam and specify output filename
$0 --bam-filter '-q 30 -f 0x2' --out test_qa.txt test.bam

## exclude chromosome chrM
$0 --out test_qa.txt --chrs 'chrM' test.bam

## run the ENCODE ATAC-seq pipeline command for single-end reads
$0 --out test_qa.txt --chrs 'chrM' --bam-filter '-F 1804' test.bam

## run the ENCODE ATAC-seq pipeline command for paired-end reads
$0 --out test_qa.txt --chrs 'chrM' --bam-filter '-F 1804 -f 2' test.bam

EOF

}

function clean_up
{
	if [[ $tmpFiles ]]; then
		rm -rf ${tmpFiles[@]}
	fi
}

function fail
{
	clean_up;
	exit $1;
}

# construct a associative array
function set_skipped_steps
{
	#declare -A skippedSteps
	skipped=($(str_split "," $1))
	i=0;
	for s in ${skipped[@]}
	do
		if [[ $s =~ ^[0-9]+$ ]]; then
			j=$(( $s - 1 ));
			skipped[$i]=${steps[$j]}; # convert number to step names
		fi
		i=$(( i + 1 ));
	done
	echo "${skipped[@]}"
}

# get read statistics from a bam, including duplication,
# chromosome-specific reads
function bam_stat
{
	bam=$1;
	filter=$2;
	chrs=$3;
	printHead=$4;

	tempBed=tmp.$$.b2b.bed
	if [[ $paired ]]; then
	#	tmpBam=$(sort_bam $bam "T")
	#	mv $tmpBam $bam
		#samtools view $bamFilter -u $bam | bamToBed \
		samtools view -u $filter $bam | bamToBed \
		-bedpe | gawk 'BEGIN{OFS="\t"}$1==$4{print $1,$2,$6,$9 "/" $10}' \
		>$tempBed
	else # single-end
		#samtools view $bamFilter -u $bam | bamToBed \
		samtools view -u $filter $bam | bamToBed \
		 | gawk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' \
		>$tempBed
	fi
	# print out header if true
	if [[ $printHead ]]; then
		header="#total_read_pairs\tdistinct_read_pairs\tnondup_read_pairs\ttwice_read_pairs\tNRF\tPBC1\tPBC2"
		for c in $(str_split "," $chrs)
		do
			header+="\t$c";
			header+="\t${c}Uniq";
		done
		echo -e "$header"
	fi
	# calculate now
	sort $tempBed | uniq -c | perl -e '
	$chrs=shift; # chromosomes to check
	$mt=$m0=$m1=$m2=0;
	while(<>) {
		@fields=split /\s+/;
		$cnt=$fields[1];
		$unique=$cnt==1?1:0;
		$mt+=$cnt; # total reads
		$m0++; # distinct reads
		$m1++ if $unique; # unique reds
		$m2++ if $fields[1]==2; # twice reads
		if($chrs=~/$fields[2](,|$)/) { # matching chromosomes
			$chr=$fields[2];
			${$chr}+=$cnt; # total for this chrom
			if($unique) {
				$var="$chr"."Uniq"; # unique ones for this chrom
				${$var}++; 
			} 
		} 
	}
	# output the summary
	$nrf=sprintf("%.6f", $m0/$mt);
	$pbc1=$m0>0? sprintf("%.6f",$m1/$m0):"NA";
	$pbc2=$m2>0? sprintf("%.6f",$m1/$m2):"NA";
	$chrInfo="";
	foreach (split /,/, $chrs)
	{
		$varT=$_;
		$varU=$_."Uniq";
		$chrInfo.="\t".(${$varT}||0)."\t".(${$varU}||0);
	}
	$chrInfo =~ s/^\t//;
	print join("\t", $mt, $m0, $m1, $m2, 
	$nrf, $pbc1, $pbc2, $chrInfo), "\n";
	' $chrs
	rm $tempBed;
}

# the order here MUST be the same as the those in usage description
steps=(markdup flagstat libcompex)

skips=""
bamFilter=""
outFile=/dev/stdout
exclChrs=""
posArgs=()


while [[ $# -gt 0 ]];
do
	k=$1; shift;
	case $k in
		--skip)
			skips=$1;
			shift;
			;;
		--bam-filter)
			bamFilter=$1;
			shift;
			;;
		--chrs)
			exclChrs=$1;
			shift;
			;;
		--out)
			outFile=$1;
			shift;
			;;
		*)
			posArgs+=($k)
			;;
	esac
done

if [[ ${#posArgs[@]} -lt 1 ]]; then
	usage;
	clean_up;
	exit 1;
fi

set -- ${posArgs[@]}
depends=(samtools picard)

for e in ${depends[@]}
do
	if [[ ! $(command -v $e) ]]; then
		msg "Command '$e' can't be found"
		fail 2
	fi
done

# get skipped steps in an associative array
declare -A skippedSteps
tmp=($(set_skipped_steps $skips))
#echo "input: $skips >> ${tmp[@]}"
for s in ${tmp[@]}
do
#	echo $s
	skippedSteps[$s]=1; # add to an associate array
done
#echo "Skipped steps are: '${!skippedSteps[@]}'"

# preprocess the string containing excluded chromosomes
exclChrs=$(echo $exclChrs | sed -e 's/ //g') # remove spaces
if [[ $exclChrs ]]; then
	exclChrs+=',';
	msg "exluded chromosome string: $exclChrs"
fi

bam=$1;
# make a copy of the input bam file
#bam=tmp.$$.copy.bam
#msg "Making a copy of input bam file $inBam to $bam"
paired=$(is_bam_paired $bam)

stepI=1;

msg "Step $stepI: mark read duplicates on $bam"
step=markdup
echo ">>Step $step" >$outFile
if [[ ${skippedSteps[$step]} ]]; then
	msg "Step $step is skipped"
else
	markedBam=${bam/%.bam/.dup_marked.bam}
	metric=${bam/%.bam/.dup_metrics.txt}
	if [[ ! $(is_bam_sorted $bam) ]]; then
		tmp=$(sort_bam $bam)
		mv $tmp $bam
	fi
	picard MarkDuplicates \
	I=$bam ASO=coordinate M=$metric O=$markedBam \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 TMP_DIR=. \
	TAGGING_POLICY=All
	mv $markedBam $bam; # change original bam to marked one
	msg "$bam is dup-marked now"
	cat $metric >>$outFile
	rm $metric
fi

if [[ $bamFilter ]]; then
	filteredBam=${bam/%.bam/.filtered.bam}
	msg "Filtering bam file $bam into $filteredBam"
	samtools view $bamFilter -b -o $filteredBam $bam
	bam=$filteredBam
	echo "# Bam file is filtered with '$bamFilter'" >>$outFile
fi

stepI=$(( stepI + 1 ))
msg "Step $stepI: run flagstat on $bam"

step=flagstat
echo ">>Step $step" >>$outFile
if [[ ${skippedSteps[$step]} ]]; then
	msg "Step $step is skipped"
else
	if [[ ! $(is_bam_sorted $bam) ]]; then
		msg "Sorting bam file $bam"
		tmp=$(sort_bam $bam)
		mv $tmp $bam
	fi
	samtools index $bam
	samtools flagstat $bam >>$outFile
fi

stepI=$(( stepI + 1 ))
msg "Step $stepI: compute library complexity on $bam"
step=libcomplex
echo ">>Step $step" >>$outFile
if [[ ${skippedSteps[$step]} ]]; then
	msg "Step $step is skipped"
else
	# sort bam if in paired-end
	if [[ $paired ]]; then
		tmpBam=$(sort_bam $bam "T")
		mv $tmpBam $bam
	fi
	bam_stat $bam "" $exclChrs "T" >>$outFile;
	mapqCut=30
	echo "# MAPQ > $mapqCut" >>$outFile;
	bam_stat $bam "-q $mapqCut" $exclChrs >>$outFile;
fi

clean_up;

msg "Job done at `date`";

exit 0;

	echo -e	"#total_read_pairs\tdistinct_read_pairs\tnondup_read_pairs\ttwice_read_pairs\tNRF\tPBC1\tPBC2" >>$outFile
	echo 
	gawk 'BEGIN{mt=0;m0=0;m1=0;m2=0}
		($2==1){m1=m1+$1}
		($2==2){m2=m2+$1}
		{mt=mt+$1*$2}
		{m0=m0+$1}
		END{
		pbc1=m0>0? sprintf("%f", m1/m0):"NA";
		pbc2=m2>0? sprintf("%f", m1/m2):"NA";
		printf "%d\t%d\t%d\t%d\t%f\t%s\t%s\n",mt,m0,m1,m2,m0/mt,pbc1,pbc2}' \
	$tempFile >>$outFile
	tmpFiles+=($tempFile)
fi

clean_up;

msg "Job done at `date`";

exit 0;

