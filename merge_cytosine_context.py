#!/usr/bin/env python

import sys;
import re;
import argparse as ap;
from datetime import datetime as dt;

## functions
def warn(txt):
    '''
    print warning message
    '''
    curTime=dt.now().strftime("%Y-%m-%d %H:%M:%S");
    print("[{0}] {1}".format(curTime, txt), file=sys.stderr);

def open_file(f):
    '''
    open a file
    '''
    openMethod=open;
    if re.search("\\.gz$", f, re.IGNORECASE):
        openMethod=gzip.open;
    h=openMethod(f,"rt");
    hObj={ "handle": h,
            "stacked": [],
            "eof": False};
    return(hObj);

def close_file(h):
    h['handle'].close();

def next_line(h):
    '''
    read next line from the file handle
    '''
    if len(h['stacked']):
        return(h['stacked'].pop(0));
    return(h['handle'].readline());

def output_record():
    '''
    output the last stored record
    '''
    global lastRecord;
    if lastRecord is not None:
        print(sep.join(lastRecord));
        lastRecord=None;

def in_same_context(rec1, rec2):
    '''
    check whether two recrods are from the same CpG/CHG context
    '''
    if rec1[strandCol] == rec2[strandCol]: # same strand
        return False;
    if rec1[chrCol] != rec2[chrCol]: # different chromosomes
        return False;
    if abs(int(rec1[posCol])-int(rec2[posCol])) != contextWidth:
        return False;
    return True;

desc='''

This program merges the cytosines from the + and - strands in the same
CpG or CHG context into one, so two lines in an input file will be
output as one line.

The program can also handle the cases where the input file contains
both merged (one line per CpG/CHG) and unmerged (two lines) cytosines.

Options:

'''

authorInfo='''
Author: Zhenguo Zhang
Email: zzhang@zymoresearch.com
''';

op=ap.ArgumentParser(
        description=desc,
        formatter_class=ap.RawTextHelpFormatter,
        epilog=authorInfo
        );

op.add_argument("--type",
        help="provide the cytosine context, choosing from CpG and CHG",
        dest="cType",
        choices=["CpG","CHG"],
        required=True
        );

op.add_argument("--sep",
        help="the field separator for input file [%(default)s]",
        dest="sep",
        default=","
        );

op.add_argument("inFile",
        help="input file"
        );

args=op.parse_args();

cType=args.cType;
sep=args.sep;

## global variables
chrCol=0;
posCol=1;
strandCol=2;
if cType == "CpG":
    contextWidth=1;
elif cType == "CHG":
    contextWidth=2;
else:
    warn("Cytosine type '{cType}' isn't supported");

inFh=open_file(args.inFile);

lastRecord=None;
counter=0;
while True:
    line=next_line(inFh);
    if line == '':
        output_record();
        break;
    counter+=1;
    if counter % 10000 == 0:
        warn(f"#>> Processing Line {counter}");
    fields=line.strip().split(sep);
    if lastRecord is None:
        lastRecord=fields;
        continue;
    # other check whether in the same context
    if fields[strandCol] == "-":
        if in_same_context(lastRecord,fields):
            output_record();
        else: # not in the same context
            output_record();
            lastRecord=fields;
    else: # + strand, new context
        output_record();
        lastRecord=fields;

close_file(inFh);

warn("Job done\n");

sys.exit(0);

