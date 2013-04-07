#!/usr/bin/perl
#group_extract.pl
##################################
# Written by Joe Edwards
# Sundaresan Lab
# Department of Plant Biology
# University of California - Davis
# 2012
##################################
use lib "/home/ledwards/Perl_Scripts/";
use strict; use warnings;
use Getopt::Std;

#==============================================================================
## Get the options from command line
my %opts;
getopts('g:f:r:', \%opts);


#==============================================================================
## Usage statement
my $usage = "
group_extract.pl
Last Updated: 5/6/2013

Written by Joe Edwards
PhD Student
Sundaresan Lab
Department of Plant Biology
University of California - Davis

Usage
group_extract.pl -f {fastq1} -r {fastq2}

Options
	-f:	Forward fastq file
	-r:	Reverse fastq file
	-g:	group file
";


#==============================================================================
## Make sure everything on the command line is kosher and open files
die "[ERROR] Please provide a forward fastq file$usage" unless exists $opts{f};
die "[ERROR] Please provide a reverse fastq file$usage" unless exists $opts{r};
die "[ERROR] Please provide a group file$usage" unless exists $opts{g};

open F1, "$opts{f}" or die "Cannot open forward fastq file $opts{f}\n";
open F2, "$opts{r}" or die "Cannot open reverse fastq file $opts{r}\n";
open G1, "$opts{g}" or die "Cannot open group file $opts{g}\n";

my ($f1out) = $opts{f} =~ /(.*)\.fq/;
my ($f2out) = $opts{r} =~ /(.*)\.fq/;

$f1out .= ".grouped.fq";
$f2out .= ".grouped.fq";

open F1OUT, ">$f1out" or die "Cannot write to file $f1out\n";
print "Writing to file $f1out\n";
open F2OUT, ">$f2out" or die "Cannot write to file $f2out\n";
print "Writing to file $f2out\n";

#==============================================================================
while (my $line = <G1>) {
	chomp $line;
	my ($name, $group) = split(/\t/, $line);
	FQread($name);
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Subroutines Follow
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub FQread {
	my ($name) = @_;
	while (my $header1 = <F1>) {
		my $header2 = <F2>;
		chomp $header1;
		chomp $header2;

		my $seq1 = <F1>; chomp $seq1;
		my $seq2 = <F2>; chomp $seq2;
		my $skip1 = <F1>; my $skip2 = <F2>;
		my $qual1 = <F1>; chomp $qual1;
		my $qual2 = <F2>; chomp $qual2;

		if ($header1 =~ m/$name/) {
			print F1OUT "$header1\n$seq1\n+\n$qual1\n";
			print F2OUT "$header2\n$seq2\n+\n$qual2\n";
			last;
		}
	}
}










