#!/usr/bin/env perl                                                                                                                                                                                                                    

###### THIS SCRIPT OPENS THE STATS FILE FROM THE GETCONTIGS SCRIPT
###### IT FINDS THE EXONS FOR EACH GENE AND STITCHES THEM TOGETHER
###### END RESULT IS ONE EXON FILE FOR EACH GENE WITH cDNA 
###### THESE OUTPUT FILES *.Exonerate_Exons.fasta  GO DIRECTLY TO ALIGNMENT


##### LIST ALL THE GENES
#`ls -l *.results.ed.fasta >result.files`;
#open FHx, "<result.files";

system "ls -l *.stats.csv >files";
open FH, "<files";
while (<FH>) {
	if (/(\S+).results.ed.stats.csv/) {
		my $gene=$1;
		open FH1, "<$gene.results.ed.stats.csv";
		open OUT1, ">$gene.Exonerate_Exons.fasta";
		while (<FH1>) {	
			#### GET THE WHOLE CONTIG IF HAS THE WHOLE GENE
			if (/^(\S+?),.*TRUE,\S+,,,,,,,,,,,,(.*)$/) {
				my $lib=$1;
				my $contig=$2;
				#print OUT1 "\n>$lib_$contig\n";
				#print "$lib\tTRUE and contig is $contig\n";
				my $seqflag=0;
				open FASTA, "<$gene.results.ed.fasta";
				while (<FASTA>) {
					$line=$_;
					chomp $line;
					if ($line =~ /^>/) {  $seqflag=0; }
					if ($seqflag == 1 ) { print OUT1 $line;}  # print $line;}
					if ($line =~ m/>$lib.*$contig/g) {
						$seqflag =1;
						print OUT1 ">$lib\n";
						#print OUT1 ">TRUE$lib\_$contig\n";
						#print  ">$lib\_$contig\n";
					}
				}
				print OUT1 "\n";
			}
			#### DOES NOT HAVE THE WHOLE CONTIG 
			else {
				my $line=$_;
				my @array = split (/,/, $line);
				#### FIND OUT IF ONLY NEED ONE CONTIG
				my $numcontigs = $array[11];
				if ($numcontigs == 1) {
					my $lib = $array[0];
					my $start=$array[14];
					my $end = $array[15];
					my $length=$array[2];
					#print "\n$gene\t$lib START $start END $end LENGTH $length\n";
					my $contig = $array[16];
					my $seqflag=0;
					my $sequence=();
					open FASTA, "<$gene.results.ed.fasta";
					while (<FASTA>) {
						$line=$_;
						chomp $line;
						if ($line =~ /^>/) { $seqflag=0; }
						if ($seqflag == 1 ) { $sequence = $sequence.$line;}
						if ($line =~ m/>$lib.*$contig$start$end/g) {
							$seqflag =1;  #print OUT1; 
						#	print OUT1 ">ONECONTIG$lib\_$contig\_$start\_$end\n";
							print OUT1 ">$lib\n";
						#print;  #"\n>$lib\_$contig\n";
						}
					}
					###### PRINTS OUT NNNs AT THE BEGINNING OF THE GENE IF THE CONTIG DOES NOT START AT THE BEGINNING
					if ($start != 0) { for (0..$start) { print OUT1 "NNN"; } }
					print OUT1 "$sequence"; #print "$sequence";
					##### PRINTS OUT NNNs AT THE END OF THE SEQUENCE IF CONTIG DOES NOT COVER THE FULL LENGTH
					if ($end != $length) { for ($end .. ($length-1)) { print OUT1 "NNN";  } } 
					print OUT1 "\n";
				}
				######### ELSE IF MORE THAN ONE CONTIG STITCH THEM TOGETHER WITH NNNS IN THE MIDDLE 
				else {
					my $lib = $array[0];
					my $length = $array[2];
					$count=0;
					$contigstart=14;
					my $gapstart=0;
					for (1..$numcontigs) {
						my $contignumber=$_;
						print "Num contigs $numcontigs position $contignumber\n";
						my $start=$array[14+$count];
						my $end=$array[15+$count];			
						my $pos=14 + (($numcontigs*2) + ($contignumber-1));
						print "Position $pos\n";
						my $contig=$array[$pos];
						my $last=0;
						my $gapend=0;
						#my $gapstart=0;
						my $sequence='';
						print "CONTIG start $start end $end contig $contig\n";
						##### FIRST CONTIG
						if ( $contignumber == 1) { 
							print "$lib FIRST\t$start\t$end\n";
							print "contig $gene\t$contig\n";
							$gapstart=$end+1;
							### REMOVE
							#print OUT1 "\n";
							### REMOVE
							open FASTA, "<$gene.results.ed.fasta";
							while (<FASTA>) {
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line; }
								if ($line =~ m/>$lib.*$contig$start$end/g) {
									$seqflag =1;  #print OUT1; 
									print OUT1 ">$lib\n";
									#print OUT1 ">FIRST$lib\_$contig\_$start\_$end\n";
									#print;  #"\n>$lib\_$contig\n";
								}
							}
							#### AGAIN PRINT OUT NNNs IF NECESSARY AT THE BEGINNING
							if ($start != 0) { 
								for (0..$start) { print OUT1 "NNN";} 
							}
							print OUT1 "$sequence";
							$nextstart=$end;
							close FASTA;
						}
						$sequence='';
						### LAST CONTIG
						if ($contignumber == $numcontigs)  { 
							print "$lib LAST\t$start\t$end\tgapstart $gapstart to $start - 1 \n";
							print "contig $gene\t$contig\n";
							#print OUT1 "GAPSTART\n";
							for ($gapstart ..($start-1)) { print OUT1 "NNN"; i} 
							#print OUT1 "LAST SEQ\n";
							print "MATCH  $lib.$contig\_$start\_$end\n";
							open FASTA, "<$gene.results.ed.fasta";
							while (<FASTA>) {
								#print;
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line;}
								if ($line =~ m/>$lib.*$contig$start$end/) {
									$seqflag =1; 
									#print OUT1 ">LAST $lib\_$contig\_$start\_$end\n";
									#print;  #"\n>$lib\_$contig\n";
								}
							}	
							print OUT1 "$sequence";
							#########  PRINT OUT NNNs AT THE END IF SEQUENCE DOES NOT COVER FULL LENGTH
							if ($end != $length) { 
								for ($end .. ($length-1)) { print OUT1 "NNN";  } 
							}
							print OUT1 "\n";
							close FASTA;
						}
						$sequence='';
						#### MIDDLE CONTIGS
						if ($contignumber > 1 && $contignumber != $numcontigs) {
							#print "\nMIDDLE CONTIG\n";
							#print "$lib MIDDLE\t$start\t$end\tgapstart $gapstart to $start - 1 \n";
							#print "$lib\t$contig\n";
							#print OUT1 "GAPSTART\n";
							##### PRINT NNNs IN BETWEEN CONTIGS
							for ($gapstart ..($start-1)) { print OUT1 "NNN";} 
							#print OUT1 "LAST SEQ\n";
							#print "MATCH  $lib.$contig\_$start\_$end\n";
							open FASTA, "<$gene.results.ed.fasta";
							$gapstart = $end + 1;
							while (<FASTA>) {
								#print;
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line;}
								if ($line =~ m/>$lib.*$contig$start$end/) {
									$seqflag =1; 
									#print OUT1 ">LAST $lib\_$contig\_$start\_$end\n";
									#print;  #"\n>$lib\_$contig\n";
								}
							}	
							print OUT1 "$sequence";
						}
						$count=$count+2;
					}
				}
			}
		}
	}
}

