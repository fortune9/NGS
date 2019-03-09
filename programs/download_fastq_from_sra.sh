#!/bin/bash

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

options="--gzip --read-filter pass --skip-technical --split-3 --clip"

if [[ $2 != '' ]]; then
	options=${@:2}
fi

cmd="fastq-dump $options $1"

echo "# Running command '$cmd' at `date`" >&2

$cmd

echo "# Downloading $1 finished at `date`" >&2

exit 0;

