#!/bin/bash

set -eu

function usage()
{
	cat <<EOF
Usage: $0 <num/frac> <outfile-base> <fq1-file> [<fq2-file>]

This program samples the reads from input files using seqtk tools.
It accepts either one or two fastq files. When two files are fed,
the sampling is down in paired-end mode.

The number of sampled reads/pairs is specified by 1st parameter
num/frac, i.e, either a number greater than 1 or a fraction smaller
than 1.

The paramter <outfile-base> specifies basename for output files. In
single end mode, fq.gz is appended, while in paired-end mode, R1.fq.gz
and R2.fq.gz are appended for the two output files, respectively.

E.g.: 
$0 0.1 out fq1.gz # single-end mode, yield out.F0_1.fq.gz

$0 0.1 out fq1.gz fq2.gz # paired-end mode, yield out.F0_1.R1.fq.gz and out.F0_1.R2.fq.gz
EOF
}

function msg
{
	echo "$*" >&2
}

if [[ $# -lt 3 ]]; then
	usage
	exit 1;
fi

depends=(seqtk)

for exe in ${depends[@]}
do
	if [[ ! $(command -v $exe) ]]; then
		echo "command '$exe' can't be found" >&2
		exit 2;
	fi
done

frac=$1;
outBase=$2;
fq1=$3;
fq2="";
twopass="";
seed=$$;

if [[ ${4:-} ]]; then
	fq2=$4;
fi

if [[ $( echo "$frac > 1" | bc ) > 0 ]]; then
	# this is a number
	echo "$frac reads will be sampled from each file with seed $seed" >&2
	twopass="-2"
else
	echo "$frac fraction of reads will be sampled from each file with seed $seed" >&2
fi

msg "Job is started at `date`"

if [[ $fq2 ]]; then
	msg paired-end mode sampling.
	out1=$outBase.R1.fq.gz
	out2=$outBase.R2.fq.gz
	seqtk sample $twopass -s $seed $fq1 $frac | gzip -c >$out1
	seqtk sample $twopass -s $seed $fq2 $frac | gzip -c >$out2
else
	out=$outBase.fq.gz
	seqtk sample $twopass -s $seed $fq1 $frac | gzip -c >$out
	msg "File $out is generated"
fi

msg "Job is done at `date`"

exit 0;

