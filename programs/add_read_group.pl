#!/usr/bin/env perl
use strict;
use Getopt::Long;

my $headFile;
my $rgFromName;

GetOptions(
	"head=s"  => \$headFile,
	"rg-name!" => \$rgFromName
);

my $inFile = shift or &usage();

# read into head information
my $headLines = '';
my $rgID;
open(H,"< $headFile") or die "Cannot open $headFile:$!";
while(<H>)
{
	next unless /^\@/;
	$headLines .= $_;
	next unless /^\@RG/; # RG line only for RG ID
	unless(defined $rgID)
	{
		($rgID) = /\tID:(\S+)\s/;
	}
}
close H;

open(IN, "< $inFile") or die "Cannot open $inFile:$!";
my $counter = 0;
while(<IN>)
{
	print and next if /^\@/; # output old head lines directly
	# insert user's head lines
	print $headLines;
	# now output this line, and give the control to the subroutines
	$counter++;
	if($rgFromName)
	{
		chomp;
		($rgID) = /^([a-zA-Z0-9]+)/o;
		print "$_\tRG:Z:$rgID\n";
		add_rg_from_name(*IN);
	}else
	{
		unless($rgID)
		{
			die "No RG ID was found in $headFile:$!";
		}
		chomp;
		print "$_\tRG:Z:$rgID\n";
		add_rg_from_user(*IN);
	}

}
close IN;


exit 0;

sub add_rg_from_name
{
	my $fh = shift;
	while(<$fh>)
	{
		chomp;
		($rgID) = /^([a-zA-Z0-9]+)/o;
		print "$_\tRG:Z:$rgID\n";
		warn "$counter read lines have been processed\n"
		if(++$counter % 10000 == 0);
	}
}

sub add_rg_from_user
{
	my $fh = shift;
	while(<$fh>)
	{
		chomp;
		print "$_\tRG:Z:$rgID\n";
	}
	warn "$counter read lines have been processed\n"
	if(++$counter % 10000 == 0);
}

sub usage
{
	print <<USAGE;
Usage: $0 [options] <sam-file>

This program adds read group information at the header section and for
each read. The read group information is often needed by the GATK tools.

<sam-file> can be given as '-', indicating the standard input.

Options:

--head: a file containing the read group information to be added at
the SAM file's header section, like:
	
	\@RG\\tID:group1\\tSM:sample1\\tPL:illumina\\tLB:lib1\\tPU:unit1
	\@RG\\tID:group2\\tSM:sample1\\tPL:illumina\\tLB:lib1\\tPU:unit2

from which the RG id in the first line will be parsed and applied to 
all sequence reads in the input SAM file, but also see below.

--rg-name: a switch option. If provided, the RG id will be parsed from
the name of each sequence read. This will be useful when the SAM file
contains multiple read groups and the read names contain relevant read
group information. At present, the RG id in each read name is assumed
to be the longest beginning string containing only alphabet
characters, i.e., matching the pattern '^[a-zA-Z0-9]' in Perl.

Author: Zhenguo Zhang
Created: Sun Nov  6 14:47:25 EST 2016

USAGE

	exit 1;
}

__END__

=pod

=head1 NAME

B<add_read_group.pl>

This is a program to add read group information to SAM read alignment
files.

=head1 LICENSE AND COPYRIGHT

Copyright 2016 Zhenguo Zhang.

GNU License 3.0 or newer

=cut
