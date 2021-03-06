#!/usr/bin/env perl                                                                                                                                                                                                                    

###############  THIS FILE GOES THROUGH THE ATRAM CONTIGS NAD STICHES THEM TOGETHER IF NECESSARY

my $overlap = shift;
my $aTRAMfile = shift;
my $path = shift;

########## EDIT THESE ROWS
#my $overlap = 10;  ## leave this at the number set for getcontigs.first.pl, our default is 10
#my $aTRAMfile = 'I4.trinity.out.best.ed';   ## this is the file name of the aTRAM output, minus the variable names
#my $path = '/data/anoplura/HaeleCombined/1107_Contigs/BESTFILES';  ## path to where the aTRAM contigs are
##############################################



open CONTIGINFO, ">Gene.Contig.Info.txt";
print CONTIGINFO "Gene\tLib\tNumcontigs\tUnique\tDiff\n";

system "ls -l *.stats.OVERLAP.$overlap.csv >files";
open FH, "<files";
while (<FH>) {
	if (/(\S+).results.stats.OVERLAP.$overlap.csv/) {
		my $gene=$1;
		open FH1, "<$gene.results.stats.OVERLAP.$overlap.csv";
		open OUT1, ">$gene.aTRAM.Exonerate.Round1.OVERLAP.$overlap.fasta";
		while (<FH1>) {
			#print;	
			#### GET THE WHOLE CONTIG IF HAS THE WHOLE GENE
			##### EDITS APRIL 3 ONLY GET THE OPEN READING FRAME
			#Btmac,3,313,TRUE,trinity,,,,,,,,0,0,313,3.c
#			$assembler,,,,,,,,0,$beg,$end,$contighash{$tax}
			if (/^(\S+?),.*TRUE,\S+,,,,,,,,0,\d+,\d+,CONTIGS,(.*)$/) {
				my $lib=$1;
				my $contig=$2;
				#print OUT1 ">$lib\_$gene\_Whole_Gene\n";
				#print "$lib\tTRUE and contig is $contig\n";
				my $seqflag=0;
				print "MY LIB IS $lib\n MY FILE IS $path\/$gene.$lib.$aTRAMfile.fasta\n\n";
				open FASTA, "<$path\/$gene.$lib.$aTRAMfile.fasta";
				while (<FASTA>) {
					$line=$_;
					#print;
					chomp $line;
					if ($line =~ /^>/) {  $seqflag=0; }
					if ($seqflag == 1 ) { print OUT1 $line;}  # print $line;}
					if ($line =~ m/>$contig/g) {
						$seqflag =1;
						print "MATCH $contig\n";
						print OUT1 ">$lib\_$gene\_whole_gene\n";
						#print OUT1 ">TRUE$lib\_$contig\n";
						#print  ">$lib\_$contig\n";
					}
				}
				print OUT1 "\n";
			}
			#### DOES NOT HAVE THE WHOLE CONTIG 
			elsif (/(\S+)/ && ! /Statistics/ && ! /Inputfile/ && ! /Allowing/ && ! /There\s+are/ && ! /Library/) {
				print "LINE 54 not full contig\n";
				my $line=$1;
				#print "$line\n";
				my @array = split (/,/, $line);
				#### FIND OUT IF ONLY NEED ONE CONTIG
				#Hbarb,2,313,FALSE,NA,trinity,1,2,1,trinity,310,0,3,313,1c
				my $numcontigs = $array[8];
				print "number of contigs is $numcontigs\n";
				if ($numcontigs == 1) {
					my $lib = $array[0];
					print "line 61  lib = $lib\n num contigs == 1\n";
					my $start=$array[12];
					my $end = $array[13];
					my $length=$array[2];
					#print "\n$gene\t$lib ONE CONTIG $numcontigs START $start END $end LENGTH $length\n";
					my $contig = $array[15];
					chomp $contig;
					#print "$contig\n";
					my $seqflag=0;
					my $sequence=();
				        print "MY FILE IS $path\/$gene.$lib.$aTRAMfile.fasta\n\n";
                               		open FASTA, "<$path\/$gene.$lib.$aTRAMfile.fasta";
					#open FASTA, "<$gene.results.ed.fasta";
					while (<FASTA>) {
						$line=$_;
						#print;
						chomp $line;
						if ($line =~ /^>/) { $seqflag=0; }
						if ($seqflag == 1 ) { $sequence = $sequence.$line;}
						if ($line =~ m/^>$contig/) {
#						if ($line =~ m/^>$lib/) {
#						if ($line =~ m/$contig/) {
							$seqflag =1; 
							print "ONE CONTIG MATCH $lib  $contig\n";
							 #print OUT1; 
							print OUT1 ">$lib\_$gene\_ONE_CONTIG\n";
						#	print OUT1 ">$lib\n";
						#print;  #"\n>$lib\_$contig\n";
						}
					}
					###### PRINTS OUT NNNs AT THE BEGINNING OF THE GENE IF THE CONTIG DOES NOT START AT THE BEGINNING
					#if ($start > 0) { for (0..$start) { print OUT1 "NNN"; } }
					print OUT1 "$sequence"; #print "$sequence";
					##### PRINTS OUT NNNs AT THE END OF THE SEQUENCE IF CONTIG DOES NOT COVER THE FULL LENGTH
					#if ($end < $length) { for ($end .. ($length-1)) { print OUT1 "NNN";  } } 
					print OUT1 "\n";
				}
				######### ELSE IF MORE THAN ONE CONTIG STITCH THEM TOGETHER WITH NNNS IN THE MIDDLE 
				else {
					my $lib = $array[0];
					print "\nlibrary $lib  more than one contig line 99\n";
					my $length = $array[2];
					$count=0;
					$contigstart=12;
					my $gapstart=0;
					#print "LINE 102\n\n should be the whole line\n $line\n";
					$line =~ s/^\S+CONTIGS,(\S+)/$1/g;
					#print "LINE 103 this should be the string of contigs\n$line\n";
					my %uniquecontighash=();
					my $countunique=0;
					my @newcontigarray= split(/,/, $line);
					for my $cont (@newcontigarray) { 
						if (! exists $uniquecontighash{$cont} ) {
								$uniquecontighash{$cont}=1;
								$countunique++;
						}
					}
					my $diff = $numcontigs - $countunique;
					#print CONTIGINFO " diff = $numcontigs - $coununique\n";
					print CONTIGINFO "$gene\t$lib\t$numcontigs\t$countunique\t$diff\n";
					for (0..($numcontigs-1)) {
						my $contignumber=$_;
						#$libcontnumberhash{$lib}=$numcontigs;
						#my %libcontnumberunique=();
						#print "line 113 Num contigs $numcontigs position $contignumber\n";
						my $start=$array[12+$count];
						my $end=$array[13+$count];			
						#my $pos=14 + (($numcontigs*2) + ($contignumber-1));
						#print "Position $pos\n";
						# REWRITING FOR NEW CONTIG ARRAY
						my $contig = $newcontigarray[$contignumber];
#						my $contig=$array[$pos];
#						my $contig=$array[$pos];
						#print "$start,$end,$contig\n";
						my $last=0;
						my $gapend=0;
						#my $gapstart=0;
						my $sequence='';
						#print "CONTIG start $start end $end contig $contig\n";
						##### FIRST CONTIG
						if ( $contignumber == 0) { 
							print "$lib FIRST\t$start\t$end\n";
							print "contig $contig\n";
							$gapstart=$end+1;
							### REMOVE
							#print OUT1 "\n";
							### REMOVE
				        		print "MY FILE IS $path\/$gene.$lib.$aTRAMfile.fasta\n\n";
 							open FASTA, "<$path\/$gene.$lib.$aTRAMfile.fasta";
							#open FASTA, "<$gene.results.ed.fasta";
							while (<FASTA>) {
								#print;
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line; }
								if ($line =~ m/>$contig/g) {
									$seqflag =1;  #print OUT1; 
									print OUT1 ">$lib\_$gene\_$numcontigs\_contigs\n";
									#print OUT1 ">FIRST$lib\_$contig\_$start\_$end\n";
									print  "\n>MATCH $contig\n";
									
								}
							}
							#### AGAIN PRINT OUT NNNs IF NECESSARY AT THE BEGINNING
							#if ($start != 0) { 
							#	for (0..$start) { print OUT1 "NNN";} 
							#}
							print OUT1 "$sequence";
							$nextstart=$end;
							close FASTA;
						}
						$sequence='';
						### LAST CONTIG
						if ($contignumber == ($numcontigs-1))  { 
							print "$lib LAST\t$start\t$end\tgapstart $gapstart to $start - 1 \n";
							print "contig $contig\n";
							#print OUT1 "GAPSTART\n";
							#for ($gapstart ..($start-1)) { print OUT1 "NNN"; } 
							#print OUT1 "LAST SEQ\n";
							#print OUT1 "  $lib.$contig\_$start\_$end\n";
				        		print "MY FILE IS $path\/$gene.$lib.$aTRAMfile.fasta\n\n";
							open FASTA, "<$path\/$gene.$lib.$aTRAMfile.fasta";
							#open FASTA, "<$gene.results.ed.fasta";
							while (<FASTA>) {
								#print;
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line;}
								if ($line =~ m/>$contig/) {
									$seqflag =1; 
									#print OUT1 ">LAST $lib\_$contig\_$start\_$end\n";
									print "\n$lib  MATCH $contig\n";
								}
							}	
							print OUT1 "$sequence";
							#########  PRINT OUT NNNs AT THE END IF SEQUENCE DOES NOT COVER FULL LENGTH
							#if ($end != $length) { 
						#		for ($end .. ($length-1)) { print OUT1 "NNN";  } 
						#	}
							print OUT1 "\n";
							close FASTA;
						}
						$sequence='';
						#### MIDDLE CONTIGS
						if ($contignumber > 0 && $contignumber < ($numcontigs-1)) {
							print "\nMIDDLE CONTIG\n";
							print "$lib MIDDLE\t$start\t$end\tgapstart $gapstart to $start - 1 \n";
							print "$lib\t$contig\n";
							#print OUT1 "GAPSTART\n";
							##### PRINT NNNs IN BETWEEN CONTIGS
							#for ($gapstart ..($start-1)) { print OUT1 "NNN";} 
							#print OUT1 "LAST SEQ\n";
							#print "MATCH  $lib.$contig\_$start\_$end\n";
				        		print "MY FILE IS $path\/$gene.$lib.$aTRAMfile.fasta\n\n";
							open FASTA, "<$path\/$gene.$lib.$aTRAMfile.fasta"; # or die "Cannot find file $path\/$gene.$lib.$aTRAMfile.fasta\n";
							#open FASTA, "<$gene.results.ed.fasta";
							$gapstart = $end + 1;
							while (<FASTA>) {
								#print;
								$line=$_;
								chomp $line;
								if ($line =~ /^>/) { $seqflag=0; }
								if ($seqflag == 1 ) { $sequence = $sequence . $line;}
								if ($line =~ m/>$contig/) {
									$seqflag =1; 
									#print OUT1 ">LAST $lib\_$contig\_$start\_$end\n";
									print "\n>$lib MATCH $contig\n";
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

