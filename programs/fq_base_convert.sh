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
Usage: $0 [options] <fq-file>

This program converts the bases in a fastq file (can be gzipped).

Options (default values are in []):

--from <char>: the old base to replace. Mandatory.

--to <char>: the new base to change into. Mandatory.

--out <string>: output filename. The default is writing to screen [].

Example uses:

# C to T conversion
$0 --from C --to T test.fq.gz | gzip -c >test.C2T.fq.gz

# G to A conversion
$0 --from G --to A test.fq.gz | gzip -c >test.G2A.fq.gz

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

outFile=/dev/stdout
from="";
to="";
posArgs=()


while [[ $# -gt 0 ]];
do
	k=$1; shift;
	case $k in
		--from)
			from=$1;
			shift;
			;;
		--to)
			to=$1;
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

if [[ ! "$from" ]] || [[ ! "$to" ]]; then
	msg "The options '--from' and '--to' are required"
	fail 2;
fi

#set -- ${posArgs[@]}

fqFile=${posArgs[0]}

msg "Processing $fqFile at `date`"

open="cat"

if [[ $fqFile =~ \.gz ]]; then
	open="zcat"
fi

$open $fqFile | gawk -v f=$from -v t=$to \
'{if(/^@/) {print; getline; gsub(f,t);print; } else {print;}}' \
>$outFile

clean_up;

msg "Job done at `date`";

exit 0;

