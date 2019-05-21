#!/usr/bin/env perl

use strict;

my $fqFile = shift;
my $strs=shift or &usage();
$strs =~ s/,/|/;
my $pat = join("","(",$strs,")");
$pat = qr/$pat$/io;

#print $pat;

if($fqFile =~ /\.gz$/i)
{
	open(I, "zcat $fqFile | ") or die $!;
}else
{
	open(I, "< $fqFile") or die $!;
}

while(my $read=_next_read())
{
	if($read->[1] =~ $pat) # matching
	{
		my ($name) = $read->[0] =~ /^@(\S+)/;
		print join("\t", $name, length($read->[1])), "\n";
	}
	# otherwise do nothing
}

close I;

print STDERR "Job is done\n";

exit 0;

sub _next_read
{
	return undef if eof(I);
	my @record;
	my $i = 0;
	while(<I>) # read 4 lines
	{
		chomp;
		push @record, $_;
		last if ++$i == 4;
	}

	if($record[0] !~ /^@/ or $record[2] !~ /^\+/)
	{
		warn $record[0]." is not a valid fastq record; skipping\n";
		_next_read(); # read next one
	}
	return \@record;
}

sub usage
{
	print <<USAGE;
Usage $0 <fastq-file> <string1[,string2,...]>

This program extracts the sequence reads whose 3' ends
matching any of the provided strings.

E.g.: $0 test.fq.gz CG

Note: the input fastq file can be gzipped and the matching
is case-insensitive.

USAGE
	exit 1;
}

