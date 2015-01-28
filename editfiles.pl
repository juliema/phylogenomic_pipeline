#!usr/bin/env perl

use strict;
use warnings;

`ls -l *.sorted >files`;

#my %contighash=();
open FH, "<files";
while (<FH>) {
	if (/(\S+).csv.sorted/) {
		my $file=$1;
#		my $num=0;
		open OUT, ">$file.ed.sorted.csv";
		open FH1, "<$file.csv.sorted";
		while (<FH1>) {
			if (/(\S+),([A-Za-z]+),(\S+)$/) {
				my $first=$1;
				my $assembler=$2;
				my $contig=$3;
#				$contighash{$contig}=$num;
				$contig =~ s/[-,@+\[\]:!]//g;
				print OUT "$first,$assembler,$contig\n";
				#$num++;
			}
		}
		open FH2, "<$file.fasta";
		open OUT2, ">$file.ed.fasta";
		#$num=0;
		while (<FH2>) {
			if (/^>(\S+),([A-Za-z]+),(\S+)$/) {	
				my $lib = $1;
				my $assembler=$2;
				my $contig=$3;
				$contig =~ s/[-,@+\[\]:!]//g;
				print OUT2 ">$lib,$assembler,$contig\n";
				#$num++;
			}
	                elsif (/^\S+/ && ! /Command|Hostname|completed/) { print OUT2; }
		}
	}
}
			
