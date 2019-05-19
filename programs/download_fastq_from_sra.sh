#!/bin/bash

set -e

if [[ $# -lt 1 ]]; then
	cat << EOF
Usage: $0 <SRR run-id> [<other fastq-dump options>]

This program is a wrapper of fastq-dump in SRA-toolkit, with the
following options are added by default:

--gzip
--read-filter pass
--skip-technical
--split-3
--clip

These default options are good for fetching paired-end sequences. If
not desired, one can specify desired parameters after NCBI SRR id; in
this case, the default options would not be activated.

E.g.: 
$0 SRR390728 # default options
$0 SRR390728 -O mydir # default options are disgarded

EOF
	exit 1;
fi

options="--skip-technical --split-3"

if [[ $2 != '' ]]; then
	options=${@:2}
fi

dump='fastq-dump'

if [[ $(command -v fasterq-dump) ]]; then
	ncores=`nproc`
	dump="fasterq-dump --threads $ncores"
else
	options+=" --gzip --read-filter pass --clip"
fi

cmd="$dump $options $1"

echo "# Running command '$cmd' at `date`" >&2

$cmd

if [[ "$dump" != "fastq-dump" ]]; then
	gzip -S .gz ${1}*.fastq; # need manual compression
fi
	
echo "# Downloading $1 finished at `date`" >&2

exit 0;

