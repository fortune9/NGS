#!/usr/bin/env python3

import sys;
import os;
import argparse as ap;
import gzip;
import re;
from datetime import datetime;

self=os.path.basename(sys.argv[0]);

## functions
def report_progress():
    '''
    report how many read names have been processed
    '''
    if counter % 10000 == 0:
        warn(f"{counter} read names have been processed");

def warn(txt):
    print(f"[{self}] "+txt, file=sys.stderr);

def id_from_file():
    '''
    read read ids from a file
    return as a generator
    '''
    while True:
        id=nameFh.readline();
        if id == '': # no more ids
            break;
        yield id.strip().split()[0]; # get the first field only

def target_ids():
    '''
    get the names of target reads as a generator.
    '''
    if namesArr == None: # names from a file
        return(id_from_file());
    else:
        return (x for x in namesArr);

def input_reads():
    '''
    get the next read record from fastq file,
    return as a generator
    '''
    while True:
        line=i.readline();
        if line == '': # no more lines
            break;
        m=re.search("^@(\S+)", line);
        if m:
            id=m.group(1);
            record=[id,line];
            # also add the next 3 lines
            record.extend([i.readline(),i.readline(),i.readline()]);
            yield record;
        else:
            print("Found line '{0}' don't match read id line format".
                    format(line), file=sys.stderr);
            continue;

desc='''
This program extracts sequence reads from a fastq file.
Depending on the input, the program acts in different modes:

(1) a number range: the read records in the given range are extracted.
For example, extracting records from 201'th to 500'th. (Not
implemented)

(2) a number or fraction: a random set of reads matching the specified
size are extracted. (Not implemented)

(3) a list of read names: the reads specified by the names are
extracted.

Note:
The input fastq file need be sorted beforehand as well as the names
provided to the option --names.

Default optional values are in [].

''';

authorInfo='''
Author:  Zhenguo Zhang
Email: zhangz.sci@gmail.com
''';

op=ap.ArgumentParser(
        description=desc,
        formatter_class=ap.RawTextHelpFormatter,
        epilog=authorInfo
        );

op.add_argument("infile",
        help="the input fastq file, can be gzipped");

## mandatory options


## auxilary options
op.add_argument("--range",
        help="the range of recoreds to be extracted, given by 2 numbers separated by comma [%(default)s]",
        type=str,
        action="store",
        default="1,10");

op.add_argument("--random",
        help="an integer or a fraction, then this number/fraction of reads are randomly extracted",
        type=float);

op.add_argument("--seed",
        help="the seed used for generating random numbers when extracting a random set of reads [current pid]",
        type=int
        )

op.add_argument("--names",
        help="the read names separated by comma or a filename with them (one per line)",
        type=str,
        action="store"
        );

op.add_argument("--outfile","-o",
        help="output filename [stdout]",
        dest="outFile",
        default=sys.stdout,
        metavar="output");

args=op.parse_args();

inFile=args.infile;

## prepare the read names to extract
namesArr=None;
nameFh=None;

warn("Job is started at {0}".format(datetime.now()));

if os.path.isfile(args.names):
    nameFh=open(args.names,"rt");
else:
    namesArr=args.names.split(',');

openMethod=open;

if re.search("\\.gz$", inFile, re.IGNORECASE):
    openMethod=gzip.open;

i=openMethod(inFile, "rt");
o=open(args.outFile,"w") if args.outFile != sys.stdout else args.outFile;

ids=target_ids(); # a generator

counter=0;
skipped=0;
for id in ids:
    #print("Searching for {0}".format(id), file=sys.stderr);
    counter+=1;
    recorded=False;
    for fqRead in input_reads():
        if id == fqRead[0]: # the id is found
            o.write("".join(fqRead[1:])); recorded=True;
            break; # get a new id and new read
        elif id < fqRead[0]: # this id has no read
            #print("[{0}] The read '{1}' has no read in input fastq".
            #        format(self, id), file=sys.stderr);
            # get new id and test
            while id < fqRead[0]:
                warn(f"Skipping id {id}"); skipped+=1;
                id=next(ids);
                counter+=1;
                report_progress();
            if id == fqRead[0]: # the id is found
                o.write("".join(fqRead[1:])); recorded=True;
                break; # get a new id and new read
            # for inner if: otherwise greater id, waiting for next read
        # for outer elif: otherwise for greater id, waiting for next read
    if not recorded:
        warn(f"Skipping id {id}"); # this can be last greater id without read
        skipped+=1;
    report_progress();

if nameFh is not None:
    nameFh.close();

i.close();
if o != sys.stdout:
    o.close();

warn("Job is done at {0}\n{1} total ids\n{2} skipped".format(
        datetime.now(), counter, skipped));

sys.exit(0);


