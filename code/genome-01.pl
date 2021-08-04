#!/usr/bin/perl

use warnings;
use strict;
use List::Util qw(sum);

#### Output the SLiM markers into genomic format.
open GENOME, $ARGV[0] || die "$!";
open MS, ">$ARGV[1].01" || die "$!";

my $nMarker=$ARGV[2];

while(<GENOME>){
	if(/#/){
		next;
	}

###write ped; pop1 is 0; pop2 is 1.
	my @g1;
	push(@g1, "1") for(0..($nMarker-1));

	my @b=split;
	$g1[$_]="0" for(@b);

### only the chromosome has more than 10% region from each parents can be considered as a hybrid. This is a weak filter cretiria for the complete ancestry chromosomes.
	my $prop= sum(@g1) / $nMarker;
	if($prop>0.1 && $prop<0.9){
		foreach my $i (0..($nMarker-1)){
			print MS " $g1[$i]";
		}
		print MS "\n";
	}
}

close GENOME; close MS;
