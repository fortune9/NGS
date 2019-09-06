#!/bin/bash

set -e

tmpFiles=();

function download_tool
{
	self=`basename $1`
	out=tmp.$$.bb.$self
	wget --quiet -O $out "$1" || exit $?;
	echo $out;
}
# load the common bash file
bashFuncUrl='https://raw.githubusercontent.com/fortune9/programs/master/bash_functions.sh'
bashFuncFile=$(download_tool $bashFuncUrl)
source $bashFuncFile
ngsFuncUrl='https://raw.githubusercontent.com/fortune9/NGS/master/ngs_functions.sh'
ngsFuncFile=$(download_tool $ngsFuncUrl)
source $ngsFuncFile
bamInfoUrl=https://raw.githubusercontent.com/fortune9/NGS/master/programs/chrom_info_from_bam.sh
bamInfo=$(download_tool $bamInfoUrl);
chmod u+x $bamInfo;

tmpFiles+=("$bashFuncFile" $ngsFuncFile $bamInfo)

#pg=${BASH_SOURCE[0]}
#pgDir=`dirname $pg`

ucscDir='http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64'


function clean_up
{
	if [[ $tmpFiles ]]; then
		#echo "cleaning ${tmpFiles[@]}"
		rm ${tmpFiles[@]};
	fi
}

function usage
{
	cat <<EOF
Usage: $0 <bam-file> <out-file.bw>

This program reads a bam and output a bigwig file.

Example uses:

$0 test.bam  test.bw
EOF
}

if [[ $# -lt 2 ]]; then
	usage;
	clean_up;
	exit 1;
fi

bamFile=$1;
outFile=$2;

bam=$bamFile
bgFile=$(rand_str).bg

msg "step 1: convert bam to bedgraph"
if [[ ! $(is_bam_sorted $bam) ]]; then
	msg "Sorting $bamFile"
	bam=$(sort $bam)
	tmpFiles+=($bam)
fi

pc=""
if [[ $(is_bam_paired $bam) ]]; then
	msg "$bamFile is in paired-end mode"
	pc="-pc"
fi
bedtools genomecov -ibam $bam -bg $pc | sort -k1,1 -k2,2n >$bgFile
tmpFiles+=($bgFile)

msg "step 2: convert bedgraph to bigwig"
bg2bwUrl=$ucscDir/bedGraphToBigWig
bg2bw=$(download_tool $bg2bwUrl)
chmod u+x $bg2bw
# get the chromosome info from bam file
chrSizeFile=$(rand_str).chrom.sizes
./$bamInfo $bam >$chrSizeFile 
tmpFiles+=($chrSizeFile)
./$bg2bw $bgFile $chrSizeFile $outFile

clean_up;

msg "Job is done at `date`"

exit 0;

