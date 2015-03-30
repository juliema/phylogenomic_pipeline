#!/usr/bin/env perl

use strict;
use warnings;

my $input = shift;
my $output = shift;


open FH, "<$input";
open OUT, ">$output";
while (<FH>) {
	if (/^[+@]\S+.*\d+\.([12]).*/) {
		my $line=$_;
		my $num=$1;
		chomp $line;
		print OUT "$line\/$num\n";
	}
	else { print OUT; }
}


