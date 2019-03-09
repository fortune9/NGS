#!/bin/bash

if [[ $# -lt 1 ]]; then
	cat <<EOF
Usage: $0 <SRP accession>

This program downloads the metadata for an input NCBI SRP accession.
The output will be stored in a file in csv format.

E.g.: $0 SRP001599
EOF

	exit 1;
fi

if [[ ! $1 =~ ^SRP[0-9]+$ ]]; then
	echo "The input SRP accession '$1' is not in right format -- SRP######" >&2
	exit 2;
fi

echo Downloading metadata for $1

URL="http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$1"

out=${1}.sample_info.csv

wget --no-verbose -O $out "$URL"

exit 0;

