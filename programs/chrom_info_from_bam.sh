#!/bin/bash

if [[ $# -lt 1 ]]; then
	cat <<EOF
Usage: $0 <bam-file>

This program parses the header section of input bam file and output a
tab-delimited file in the format:
chr1	length1
chr2	length2
...	...

EOF
	exit 1;
fi

if [[ ! $(command -v samtools) ]]; then
	echo "Command 'samtools' can't be found"
	exit 2;
fi

samtools view -H $1 | grep '^@SQ' | tr ':' '\t' \
| cut -f 3,5

exit 0;

