#!/usr/bin/env perl
use strict;
#use BerkeleyDB;
use Getopt::Long;

$SIG{INT} = \&clean_tmp;
$SIG{TERM} = \&clean_tmp;

my $addSuf;
my $minLen;

GetOptions(
	'suffix!'	=> \$addSuf,
	'length:i'	=> \$minLen
);

$addSuf = 0 unless(defined $addSuf);
$minLen = 0 unless(defined $minLen);

my $minLen1=$minLen; # read 1
my $minLen2=$minLen; # read 2

my ($fq1,$fq2,$out) = @ARGV;

&usage() unless($out);

my $fh1 = _open_file($fq1);
my $fh2 = _open_file($fq2);
my $outFh1 = _open_file("${out}_1.fastq.gz", ">");
my $outFh2 = _open_file("${out}_2.fastq.gz", '>');
my $outFhUp = _open_file("${out}_up.fastq.gz", '>');


my $tmpFile = $ENV{'TEMPDIR'}."/$fq2.$$.tmp"; # fq2 itself may contain
# folder, so it is safer to add suffix only

my $counter = 0;

warn "# Filtering reads from $fq1 and $fq2 at ".scalar(localtime)."\n";

my $upCnt = 0;
my $pairCnt = 0;
my $removeCnt = 0;
while(1)
{
	my ($read1Ref, $name1, $len1) = _get_read($fh1, '1');
	last unless($read1Ref); # no more reads
	my ($read2Ref, $name2, $len2) = _get_read($fh2, '2');
	if($name1 eq $name2) # matched reads
	{
		if($len1 >= $minLen1 and $len2 >= $minLen2)
		{
			print $outFh1 $$read1Ref;
			print $outFh2 $$read2Ref;
			$pairCnt++;
		}else
		{
			$removeCnt++;
		}
	}else # names not match
	{
		print $outFhUp $$read1Ref, $$read2Ref;
		$upCnt++;
	}

	warn "[INFO] $counter read pairs have been processed\n"
	if(++$counter % 500000 == 0);
}

close $fh1;
close $fh2;
close $outFh1;
close $outFh2;
close $outFhUp;

warn "# The whole work is done\n";
warn sprintf("#Good pairs: %10d\n#Short pairs: %6d\n#Unpaired: %6d\n",
	$pairCnt, $removeCnt, $upCnt);

exit 0;

sub _seek_read
{
	my ($fh,$pos) = @_;

	my $read;
	seek($fh,$pos->[0],0); # set new start location
	read($fh,$read,$pos->[1]); # read the content
	return \$read;
}

sub _get_read
{
	my $fh = shift;
	my $suffix = shift;

	return () if eof $fh;
	$suffix = 0 unless($addSuf); # suppress suffix if necessary

	my $cnt = 0;
	my $read = '';
	my $name;
	my $length;
	my $hasSuffix; # indicate whether the sequence name already
	# contains suffix, if so, the original suffix would be kept and
	# will not add new one
	while($cnt++ < 4)
	{
		my $line = <$fh>;
		if($cnt == 1) # read name line
		{
			$hasSuffix = 1 if($line =~ /^\@\S+[\-\/]\d+\s/);
			($name) = $line =~ /^(\@\S+)/;
			if($suffix)
			{
				$line =~ s/^(\@\S+)/$1\/$suffix/ unless($hasSuffix); #
				# do not add suffix if already exists
			}

		}elsif($cnt == 3)
		{
			if($suffix)
			{
				$line =~ s/^(\+\S+)/$1\/$suffix/ unless($hasSuffix);
			}
		}elsif($cnt == 2)
		{
			my $tmp=$line;
			chomp $tmp;
			$length=length($tmp);
		}

		$read .= $line;
	}

	$name =~ s/[\-\/]\d+$// # remove the suffix from the name
	if($hasSuffix);
	return (\$read, $name, $length);
}

sub _open_file
{
	my $file = shift;
	my $type = shift;
	my $fh;

	if($file =~ /\.gz$/i)
	{
		if($type eq '>') # for output
		{
			open($fh, " | gzip >$file") or die "Can not open $file:$!";
		}else
		{
#			open($fh, "zcat $file |") or die "Can not open $file:$!";
			open($fh, "gzip -dc $file |") or die "Can not open $file:$!";
		}
	}else
	{
		open($fh, "$type $file") or die "Can not open $file:$!";
	}

	return $fh;
}

sub clean_tmp
{
	unlink "$tmpFile" if(-e "$tmpFile");
}

sub usage
{
	print <<USAGE;
Usage: $0 [options] <fastq1> <fastq2> <out>

This program is to filter paired reads from the two input files,
and output three files: <out>_1.fastq.gz <out>_2.fastq.gz
<out>_up.fastq.gz

Currently, it relies on the read name (the leading non-blank character string
at the \@ line of fastq record) to match a pair of reads.

Options:

--suffix: a switch option, if provided, suffix /1 and /2 will be added to
read names. Default is False, no suffixes are added.
--length: <int> a read pair with either read shorter than <int> will be removed.
Default is 0, so all read pairs are output.

Author: Zhenguo Zhang
Date: Wed May 22 17:50:43 EDT 2019

USAGE

	exit 1;

}

