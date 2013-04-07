#!/usr/bin/perl
#grouper.pl
##################################
# Written by Joe Edwards
# Sundaresan Lab
# Department of Plant Biology
# University of California - Davis
# 2012
##################################
use strict; use warnings;
use Getopt::Std;
use Data::Dumper;

#==============================================================================
## Get the options from command line
my %opts;
getopts('l:f:r:m:p:', \%opts);


#==============================================================================
## Usage statement
my $usage = "
grouper.pl
Last Updated: 

Written by Joe Edwards
PhD Student
Sundaresan Lab
Department of Plant Biology
University of California - Davis

Usage
grouper.pl [ OPTIONS ] -l [list] -f [forward.fq] -r [reverse.fq]

Options
	-l:	List file containing oligos associated with groups
	-f:	FastQ file containing forward reads
	-r:	FastQ file contianing reverse reads
	-m:	Number of mismatches allowed (defaults to 2)
	-p:	Prefix for output file
";


#==============================================================================
## Make sure everything on the command line is kosher and open the files
die "[ERROR] Please specify a forward fastq file\n$usage" if !exists $opts{f};
die "[ERROR] Please specify a reverse fastq file\n$usage" if !exists $opts{r};
die "[ERROR] Please specify a list file\n$usage" if !exists $opts{l};

open LIST, "$opts{l}" or die "Cannot open list file $opts{l}\n";
print "Opening $opts{l} as list file\n";
open F1, "$opts{f}" or die "Cannot open forward $opts{f}\n";
print "Opening $opts{f} as forward fastq file\n";
open F2, "$opts{r}" or die "Cannot open reverse $opts{r}\n";
print "Opening $opts{r} as reverse fastq file\n\n";

my $mm = 2; 
if (!exists $opts{m}) {
	print "No mismatches specified:  Using 2 as default\n\n";
} else {
	print "Using $opts{m} as acceptable mismatches\n\n";
	$mm = $opts{m}
}

my $output = $opts{p} . ".groups";
open OUT, ">$output" or die "Cannot write to file $output\n";
print "Writing to file $output\n";


#==============================================================================
## Read the list file
my %oligos = ();
my %foligos = ();
my %roligos = ();

while (my $line = <LIST>) {
	chomp $line;
	my ($fbar, $rbar, $sample) = split(/\t/, $line);
	$rbar = revcomp($rbar);
	$oligos{$fbar}->{$rbar} = $sample;

	$foligos{$fbar} = 1 if !exists $foligos{$fbar};
	$roligos{$rbar} = 1 if !exists $roligos{$rbar};		
}


#==============================================================================
## Read the fastqs and align the groups
my $analyzed = 0;
my $found = 0;

while (my $header1 = <F1>) {
	chomp $header1;
	my $header2 = <F2>;
	chomp $header2;

	my ($name) = $header1 =~ /.*-(.*)\s/;
	
	my $fseq = <F1>; chomp $fseq;
	my $rseq = <F2>; chomp $rseq;

	$found = $found + compare($name, $fseq, $rseq);
	
	## Skip through non important stuff
	my $skip1 = <F1>; chomp $skip1;
	$skip1 = <F1>;	chomp $skip1;
	my $skip2 = <F2>; chomp $skip2;
	$skip2 = <F2>; chomp $skip2;

	$analyzed++;
	print "$analyzed analyzed\t$found groups identified\n" if ($analyzed % 10000 == 0);
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Subroutines Follow
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Reverse Complement Generator
sub revcomp {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/ATCGatcg/TAGCtagc/;
	return $seq;
}

## How many pairwise differences between two strings
sub hamming {													
	return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

## Do the comparisons
sub compare {
	my ($name, $seq1, $seq2) = @_;
	my $forward = bestmatch($seq1, %foligos);
	my $reverse = bestmatch($seq2, %roligos);		
	
	my $value = 0;
	if ($forward eq "" || $reverse eq "") {
		return ($value);
	} else {
		my $group = $oligos{$forward}->{$reverse};	
		print OUT "$name\t$group\n";
		$value++;
		return ($value);
	}
}

## Match up the oligos
sub bestmatch {
	my ($seq, %hash) = @_;

	my $bestmatch = "";
	my $bestvalue = 12;

	foreach my $oligo (keys %hash) {
		my $value = hamming($oligo, $seq);
		if ($value < $bestvalue) {
			$bestmatch = $oligo;
			$bestvalue = $value;
		}
	}
	if ($bestvalue <= $mm) {
		$bestmatch = $bestmatch;
	} else {
		$bestmatch = "";
	}
	return $bestmatch;
}





