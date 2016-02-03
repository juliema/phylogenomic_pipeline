#!usr/bin/env perl

use strict;
use warnings;

`ls -l *.best.fasta >files`;

open FH, "<files";
while (<FH>) {
	if (/(\S+).best.fasta/) {
		my $file=$1;
		open FH2, "<$file.best.fasta";
		open OUT2, ">$file.best.ed.fasta";
		my $num=0;
		while (<FH2>) {
			if (/^>(\S+)$/) {	
				my $contig = $1;
				print "$contig\n";
				$num++;
 				$contig =~ s/[-,=@+\[\]:!]//g;
				#$contig =~ s/[\[\]]//g;
				#$contig =~ s/[,-@+\[\]:!]//g;
				print "$contig\n";
				print OUT2 ">$contig\_$num\n";
			}
			elsif (/^\S+/ && ! /Command|Hostname|completed/) { print OUT2; }
		}
	}
}

