#!/usr/local/bin/perl

use warnings;
use strict;

open IN, $ARGV[0] || die "$!";

my $tmp=<IN>; # remove markers
my $par=$ARGV[1];


while(<IN>){
	my @a=split;
### Get individual ID;
	my $file=shift(@a);
	$file=~s/_read_1//g;
	$file=~s/.fastq//g;

	open OUT, ">$ARGV[2]/$file\_$par" || die "$!";
### print individual ID to output
	print OUT "$file";
### filter ancestry informative sites with proterior probablity >0.95 homozygous for one parental population.
	@a=map{$_>0.95 ? ($par):(0)} @a;
### print ancestry informative sites
	print OUT join("\n", @a);
	close OUT;
}

close IN;
