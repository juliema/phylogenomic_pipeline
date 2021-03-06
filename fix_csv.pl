#!/usr/bin/perl

use strict;
use warnings;

my $input=shift;
my $output=shift;

open FH, "<$input";
open OUT, ">$output";

#PHUM037810

while (<FH>) {
    if (/(Library_ID\S+)/) {
        print OUT "$1\n";
    }
    if (/^PHUM\d+,(\S+.*)/) {
        print OUT "$1\n";
    }
}
