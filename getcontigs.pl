#!/usr/bin/env perl

use strict;
use warnings;

###  Contigs need to be in order they appear in the gene. so sort by the beginning of the contig.
### This script will look for all files in the folder called *.sorted, these should be exonerate results.
### It will create a file for each input file called  $input.stats.csv with all of the information about the contigs from exonerate

### rewriting script so there is only one loop for each of the multiple contigs.  so if Trinity > 1 go to same loop and if Trinity and Abyss > 1 go through the same Trinity loop.


#open STATS, ">OUTPUT.csv"; my $date=localtime(time);

`ls -l *.sorted.csv >files`;
open FILES, "<files"; 
while (<FILES>) {
	if (/(\S+).sorted.csv/) {
		my $inputfile=$1; 
		print "$inputfile\n";
		#$inputfile =~ s/.csv//g;
		my $statoutput= "$inputfile.stats.csv";
		open STATS, ">$statoutput"; my $date=localtime(time);
		print STATS "Statistics from aTRAM assemblies on $date\n";
		print STATS "Inputfile\t$inputfile.csv\n";
		# print "Statistics from aTRAM assemblies on $date\n";
		print "Inputfile\t$inputfile.sorted\n";
		my @taxarray=(); my $counttax=0; my %taxhash=();
		###############  Get list of taxa from file.
		open FH, "<$inputfile.sorted.csv";
		while (<FH>) { if (/^(\S+?),/ && ! /Library/) {  my $tax=$1; if (! exists $taxhash{$tax}) { $taxhash{$tax}=1; push @taxarray, $tax; $counttax++; } } }
		print STATS "There are $counttax Libraries in $inputfile\n";	
		print "There are $counttax Libraries in $inputfile\n";
		close FH;	
		print STATS "Library,Number of Contigs,GeneLength,Full Gene (T/F),AssemblerFullGene,Assembler1,LastContigNumber,NumberofContigs,Assembler2,LastContigNumber,TotalNumberofContigs,ContigstoKeep,Assembler,CombinedContigLength,Beginning,End,Beginning,End,ContigName\n";
		############### Now for each taxon in the file loop through and find the best contig(s).
		for my $tax (@taxarray) {
			print STATS "$tax,"; print "\n\n$tax\n\n";
			my @contigarray=(); my $numcontigs=0; my @bigarray=(); my %contighash=();
			############  count the number of contigs for each taxon, and push the line to array of arrays, and contig into an array.
			open FH1, "<$inputfile.sorted.csv";
			my $pos=0;
			while (<FH1>) {
				if (/$tax,(\S+?,\S+?,\S+?,\S+?,\S+?),(.*)$/ && ! /Library/) { 
					my $line = $1; my $contig = $2; push @contigarray, $contig; $numcontigs++; my @linearray = split(/,/, $line);
					for my $each (@linearray) { push @{$bigarray[$pos]}, $each;}
					$pos++;	
				}	
			}
			close FH1;
			print STATS "$numcontigs,$bigarray[0][0],";  ## printing the number of contigs and the length of the full gene. 
			#####  Get statistics on each contig and find the best one per taxon.
			my $flag=0; my $fullgeneT=0; my $fullgeneA=0;
			for (0..($numcontigs-1)) {
				my $cont=$_;
				my $cdslength=$bigarray[$cont][0];
				my $alcontlength=$bigarray[$cont][1];
				#### if there is a contig with the whole length of the gene take it. Preference is trinity.
				if ($cdslength == $alcontlength) { 
					if ($bigarray[$cont][4] eq 'trinity') {
						$contighash{$tax}=$contigarray[$cont];
						$fullgeneT=1;
						$fullgeneA=0;
						$flag=1;
					}
					elsif ($flag == 0) { 
						$contighash{$tax}=$contigarray[$cont]; $fullgeneA=1; $flag=1; 
					}
				}
			}
			####  if we have a full gene for either trinity or abyss print out and end.
			if ($flag == 1 ) { 
				if ($fullgeneT == 1 ) {	print "fullgene\tTrinity\n"; print STATS "TRUE,Trinity,,,,,,,,,,,,$contighash{$tax}"; } 
				if ($fullgeneA == 1 ) { print "fullgene\tAbyss\n"; print STATS "TRUE,Abyss,,,,,,,,,,,,$contighash{$tax}"; }
			}
			#### If there is not one contig that is the whole length, go to last iteration, how many contigs from each assembler?  and examine the length of the contigs.
			if ($flag == 0) {
				# These variables are needed from every loop
				my ($Asum, $Tsum, $numlastTcontig, $numlastAcontig, $Anumtokeep, $Tnumtokeep, $lastAcontig, $lastTcontig) = (0) x 8;
				my (@Tbegarray,  @Tendarray, @Abegarray, @Aendarray, @Tcontigarray, @Acontigarray, @AKeepcontigarray, @TKeepcontigarray, %AKeepcontighashBeg, %AKeepcontighashEnd, %TKeepcontighashBeg, %TKeepcontighashEnd);
				print STATS "FALSE,NA,";  ## not a full gene so FALSE in stats file and will give assembler later.
				for (0..($numcontigs-1)) {
					my $pos=$_;
					if ($bigarray[$pos][4] eq 'trinity') {  
						my $contig=$contigarray[$pos]; $contig =~ s/^(\d)\S+/$1/; 
						if ($contig == $lastTcontig) { $numlastTcontig++; }
						elsif ($contig > $lastTcontig) {  $numlastTcontig=1; $lastTcontig=$contig; }
					}
					elsif ($bigarray[$pos][4] eq 'abyss') { 
						my $contig=$contigarray[$pos]; $contig =~ s/^(\d)\S+/$1/;
						if ($contig  == $lastAcontig) { $numlastAcontig++; }
						elsif ($contig > $lastAcontig) { $numlastAcontig=1; $lastAcontig=$contig; }
					}
				}	
				print STATS "Trinity,$lastTcontig,$numlastTcontig,Abyss,$lastAcontig,$numlastAcontig,";
				print  "Trinity,Last Contig $lastTcontig,Number of contigs last iteration $numlastTcontig\n";
				print "Abyss,Last Contig $lastAcontig,Number of contigs $numlastAcontig,\n\n";
				#########################################################
				##### Trinity == 1 and Abyss == 1 
				#########################################################
				if ($numlastTcontig == 1) {
					$Tnumtokeep=1;  
					for (0..($numcontigs-1)) {
						my $pos=$_;
						if ($bigarray[$pos][4] eq 'trinity') {
							my $contig=$contigarray[$pos];
							$contig =~ s/^(\d)\S+/$1/;
							if ($contig == $lastTcontig)  {
								$Tsum=$bigarray[$pos][1];
								my $Tcontig=$contigarray[$pos]; push @TKeepcontigarray, $Tcontig;
								my $Tbeg=$bigarray[$pos][2]; push @Tbegarray, $Tbeg; $TKeepcontighashBeg{$Tcontig}=$Tbeg;
								my $Tend=$bigarray[$pos][3]; push @Tendarray, $Tend; $TKeepcontighashEnd{$Tcontig}=$Tend;
							}
						}
					}
				}
				if ($numlastAcontig == 1) { 
					$Anumtokeep=1;
					for (0..($numcontigs-1)) {
						my $pos=$_;
						if ($bigarray[$pos][4] eq 'abyss') {
							my $contig=$contigarray[$pos];
							$contig =~ s/^(\d)\S+/$1/;			
							if ($contig == $lastAcontig)  {
								$Asum=$bigarray[$pos][1];
								my $Acontig=$contigarray[$pos]; push @AKeepcontigarray, $Acontig;
								my $Abeg=$bigarray[$pos][2]; push @Abegarray, $Abeg; $AKeepcontighashBeg{$Acontig}=$Abeg;
								my $Aend=$bigarray[$pos][3]; push @Aendarray, $Aend; $AKeepcontighashEnd{$Acontig}=$Aend;
							}
						}
					}
				}
				### These variables are only needed for those that have > 1 contig
				my (@Alengtharray, @Tlengtharray); 
				#########################################################
				##### Abyss > 1  
				#########################################################
				if ($numlastAcontig > 1) {
					for (0..($numcontigs-1)) {
						my $pos=$_;
						if ($bigarray[$pos][4] eq 'abyss') {
							my $contig=$contigarray[$pos];
							$contig =~ s/^(\d)\S+/$1/;
							#$tax last abyss  contig  $lastAcontig  number of contigs $numlastAcontig\n";
							if ($contig == $lastAcontig) {
								#print "last Abyss contig $contig == $lastAcontig\n";
								my $Abeg=$bigarray[$pos][2];
								my $Aend=$bigarray[$pos][3];
								push @Abegarray, $Abeg;
								push @Aendarray, $Aend;
								push @Acontigarray, $contigarray[$pos];
								push @Alengtharray, $bigarray[$pos][1];
								print "POSITION $pos abyss length $bigarray[$pos][1]\t beginning  $Abeg\t end $Aend\n";
							}    
						}
					}
					my $pos=0;
					my $next=1;
					for (0..$numlastAcontig-1) {
						print "POSITION $pos NEXT $next\n";
						if ($next <= ($numlastAcontig-1)) {
							### Does not overlap, stitch together
							if ($Abegarray[$next] > $Aendarray[$pos]) { print "does not overlap $Abegarray[$next] > $Aendarray[$pos] LINE 179\n";
								if ($Asum == 0 ) { $Asum = $Alengtharray[$pos]+$Alengtharray[$next]; }
								elsif ($Asum > 0 ) { $Asum = $Asum +$Alengtharray[$next]; } print "abyss sum is $Asum\n";
								my $contF=$Acontigarray[$pos];
								my $contS=$Acontigarray[$next];
								if (! exists $AKeepcontighashBeg{$contF} ) {
									push @AKeepcontigarray, $contF;
									$AKeepcontighashBeg{$contF}=$Abegarray[$pos];
									$AKeepcontighashEnd{$contF}=$Aendarray[$pos];
									$Anumtokeep=$Anumtokeep+1;			
								}
								if (! exists $AKeepcontighashBeg{$contS} ) {
									push @AKeepcontigarray, $contS;
									$AKeepcontighashBeg{$contS}=$Abegarray[$next]; 
									$AKeepcontighashEnd{$contS}=$Aendarray[$next];
									$Anumtokeep=$Anumtokeep+1;			
								}
								$pos=$next;
								$next++;
							}
							### overlaps, find the longest;
							elsif ($Abegarray[$next] <= $Aendarray[$pos]) { print "overlaps $Abegarray[$next] <= $Aendarray[$pos]\n"; print "abyss sum = $Asum\n";
								if ($Asum == 0) { print "There are no contigs kept yet $Asum = 0\n";	
									my $firstsum=$Alengtharray[$pos];
									my $secondsum=$Alengtharray[$next];
									if ($firstsum >= $secondsum) { print "$firstsum >= $secondsum keep the fist contig\n";
										$Asum = $firstsum; 
										my $cont=$Acontigarray[$pos];
										push @AKeepcontigarray,  $cont;   #$Acontigarray[$pos];
										$AKeepcontighashBeg{$cont} = $Abegarray[$pos];  #print " KEEP FIRST ADDING BEGINNING $Abegarray[$pos]\t";
										$AKeepcontighashEnd{$cont} = $Aendarray[$pos];  #print " KEEP FIRST ADDING END  $Aendarray[$pos]\n";
										$Anumtokeep=$Anumtokeep+1;
									}
									elsif ($firstsum < $secondsum) { print "$firstsum < $secondsum  keep second contig\n";
										$Asum = $secondsum;
										######## KEEP SECOND CONTIG	
										my $cont=$Acontigarray[$next];
										push @AKeepcontigarray, $cont;   #$Acontigarray[$next];
										$AKeepcontighashBeg{$cont} = $Abegarray[$next]; # print " KEEP SECOND ADDING BEGINNING $Abegarray[$next]\t";
										$AKeepcontighashEnd{$cont} = $Aendarray[$next]; # print " KEEP SECOND  ADDING END $Aendarray[$next]\n";
										$Anumtokeep=$Anumtokeep+1;
										$pos=$next;
									}
									$next++;
								}
								elsif ($Asum != 0) { print "Already have contigs in memory\n";
									my $secondsum=$Alengtharray[$next];
									if ($Asum >= $secondsum) {
										#$Askipsecond = 1;
										### keep the first 
										$next++;
										print "$Asum >= $secondsum keep original contigs\n";
										#$Asum = $Asum;
										# print "Asum != 0  && Asum >= secondsum abyss sum = $Asum\n";
									}
									elsif ($Asum < $secondsum) { print "$Asum < $secondsum keep second contig\n";
										$Asum = $secondsum; 
										#print "Asum != 0  && Asum < secondsum abyss sum = $Asum\n";
										my $cont=$Acontigarray[$next];
										undef (@AKeepcontigarray);
										undef (%AKeepcontighashBeg);
										undef (%AKeepcontighashEnd);
										push @AKeepcontigarray, $Acontigarray[$next]; 
										$AKeepcontighashBeg{$cont} = $Abegarray[$next]; # print  " KEEP SECOND ADDING BEGINNING $Tbegarray[$next]\t";
										$AKeepcontighashEnd{$cont} = $Aendarray[$next]; # print " KEEP SECOND  ADDING END $Tendarray[$next]\n";
										$Anumtokeep=1;
										$pos=$next;
										$next++;
									}
								}			    		
							}
						}
					}
				}
				$Anumtokeep=scalar(@AKeepcontigarray); print "$Anumtokeep contigs to keep\n"; print "aybss sum = $Asum\n";
				for (0..$Anumtokeep-1) { $pos=$_; print "KEEPING $pos $AKeepcontigarray[$pos]\n"; } 
				#################################################################################
				#     Trinity > 1
				################################################################################
				if ($numlastTcontig > 1) {
					for (0..($numcontigs-1)) {
						my $pos=$_;
						if ($bigarray[$pos][4] eq 'trinity') {
							my $contig=$contigarray[$pos];
							$contig =~ s/^(\d)\S+/$1/;
							if ($contig == $lastTcontig) {
								#print "last Abyss contig $contig == $lastAcontig\n";
								my $Tbeg=$bigarray[$pos][2];
								my $Tend=$bigarray[$pos][3];
								push @Tbegarray, $Tbeg;
								push @Tendarray, $Tend;
								push @Tcontigarray, $contigarray[$pos];
								push @Tlengtharray, $bigarray[$pos][1];
								print "POSITION $pos trinity length $bigarray[$pos][1]\t beginning  $Tbeg\t end $Tend\n";
							}
						}
					}
					my $pos=0;
					my $next=1;
					for (0..$numlastTcontig-1) {
						print "POSITION $pos NEXT $next\n";
						if ($next <= ($numlastTcontig-1)) {
							### Does not overlap, stitch together
							if ($Tbegarray[$next] > $Tendarray[$pos]) { print "does not overlap $Tbegarray[$next] > $Tendarray[$pos] LINE 301\n";
								if ($Tsum == 0 ) { $Tsum = $Tlengtharray[$pos]+$Tlengtharray[$next]; }
								elsif ($Tsum > 0 ) { $Tsum = $Tsum +$Tlengtharray[$next]; } print "trinity sum is $Tsum\n";
								my $contF=$Tcontigarray[$pos];
								my $contS=$Tcontigarray[$next];
								if (! exists $TKeepcontighashBeg{$contF} ) {
									push @TKeepcontigarray, $contF;
									$TKeepcontighashBeg{$contF}=$Tbegarray[$pos]; 
									$TKeepcontighashEnd{$contF}=$Tendarray[$pos];
									$Tnumtokeep=$Tnumtokeep+1;			
								}
								if (! exists $TKeepcontighashBeg{$contS} ) {
									push @TKeepcontigarray, $contS;
									$TKeepcontighashBeg{$contS}=$Tbegarray[$next]; # print "KEEP BOTH ADDING BEGINNING $Abegarray[$pos]\t";
									$TKeepcontighashEnd{$contS}=$Tendarray[$next]; # print " KEEP BOTH ADDING END  $Aendarray[$next]\n\n\n";
									$Tnumtokeep=$Tnumtokeep+1;			
								}
								$pos=$next;
								$next++;
							}
							### overlaps, find the longest;
							elsif ($Tbegarray[$next] <= $Tendarray[$pos]) { print "overlaps $Tbegarray[$next] <= $Tendarray[$pos]\n"; print "trinity sum = $Tsum\n";
								if ($Tsum == 0) { print "There are no contigs kept yet $Tsum = 0\n";	
									my $firstsum=$Tlengtharray[$pos];
									my $secondsum=$Tlengtharray[$next];
									if ($firstsum >= $secondsum) { print "$firstsum >= $secondsum keep the fist contig\n";
										#$Askipfirst = 1;
										$Tsum = $firstsum; 
										#### POS STAYS THE SAME KEEPING FIRST
										#####################################################################
										### NEEED	A FLAG TO INDICATE IF WE SKIP A CONTIG AND WHICH ONE #####
										######################################################################
										my $cont=$Tcontigarray[$pos];
										push @TKeepcontigarray,  $cont;   #$Acontigarray[$pos];
										$TKeepcontighashBeg{$cont} = $Tbegarray[$pos];  #print " KEEP FIRST ADDING BEGINNING $Abegarray[$pos]\t";
										$TKeepcontighashEnd{$cont} = $Tendarray[$pos];  #print " KEEP FIRST ADDING END  $Aendarray[$pos]\n";
										$Tnumtokeep=$Tnumtokeep+1;
									}
									elsif ($firstsum < $secondsum) { print "$firstsum < $secondsum  keep second contig\n";
										$Tsum = $secondsum;
										######## KEEP SECOND CONTIG
										
										my $cont=$Tcontigarray[$next];
										push @TKeepcontigarray, $cont;   #$Acontigarray[$next];
										$TKeepcontighashBeg{$cont} = $Tbegarray[$next]; # print " KEEP SECOND ADDING BEGINNING $Abegarray[$next]\t";
										$TKeepcontighashEnd{$cont} = $Tendarray[$next]; # print " KEEP SECOND  ADDING END $Aendarray[$next]\n";
										$Tnumtokeep=$Tnumtokeep+1;
										$pos=$next;
									}
									$next++;
								}
								elsif ($Tsum != 0) { print "Already have contigs in memory\n";
									my $secondsum=$Tlengtharray[$next];
									if ($Tsum >= $secondsum) {
										#$Askipsecond = 1;
										### keep the first 
										$next++;
										print "$Tsum >= $secondsum keep original contigs\n";
										#$Asum = $Asum;
										# print "Asum != 0  && Asum >= secondsum abyss sum = $Asum\n";
									}
									elsif ($Tsum < $secondsum) { print "$Tsum < $secondsum keep second contig\n";
										$Tsum = $secondsum; 
										#print "Asum != 0  && Asum < secondsum abyss sum = $Asum\n";
										my $cont=$Tcontigarray[$next];
										undef (@TKeepcontigarray);
										undef (%TKeepcontighashBeg);
										undef (%TKeepcontighashEnd);
										push @TKeepcontigarray, $Tcontigarray[$next]; 
										$TKeepcontighashBeg{$cont} = $Tbegarray[$next]; # print  " KEEP SECOND ADDING BEGINNING $Tbegarray[$next]\t";
										$TKeepcontighashEnd{$cont} = $Tendarray[$next]; # print " KEEP SECOND  ADDING END $Tendarray[$next]\n";
										$Tnumtokeep=1;
										$pos=$next;
										$next++;
									}
								}			    		
							}
						}
					}
				}
				$Tnumtokeep=scalar(@TKeepcontigarray); print "$Tnumtokeep contigs to keep\n"; print "trinity sum = $Tsum\n";
				for (0..$Tnumtokeep-1) { $pos=$_; print "KEEPING $pos $TKeepcontigarray[$pos]\n"; } 
				##################################################################
				#          FINAL COMPARISON
				##################################################################
			    	if ($Asum > $Tsum ) { 
					print "Keep abyss contigs $Asum > $Tsum\n\n\n\n";
					print STATS "$Anumtokeep,Abyss,$Asum,";
					for my $contig  (@AKeepcontigarray) {
			    			print STATS "$AKeepcontighashBeg{$contig},$AKeepcontighashEnd{$contig},";
			    			print  "sum $Asum > $Tsum keeping Abyss $AKeepcontighashBeg{$contig},$AKeepcontighashEnd{$contig},$contig,";
					}
					for my $contig  (@AKeepcontigarray) {
			    			print STATS "$contig,";
					}
		    		}
		    		elsif ($Asum <= $Tsum ) { 
					print "Keep trinity contigs $Asum <= $Tsum\n\n\n";
					print  "number to keep $Tnumtokeep,Trinity,";
					print STATS "$Tnumtokeep,Trinity,$Tsum,";
					for my $contig  (@TKeepcontigarray) { 
			    			print STATS "$TKeepcontighashBeg{$contig},$TKeepcontighashEnd{$contig},";  	
					}
					for my $contig  (@TKeepcontigarray) { 
			    			print STATS "$contig,";  	
			    			print  "$Tsum > $Asum\n$TKeepcontighashBeg{$contig},$TKeepcontighashEnd{$contig},$contig\n"; 
					}
		    		}
			}	    
			##################################################################
			#          END!!!
			##################################################################
			print STATS "\n";
			undef(@contigarray);
			undef(@bigarray);
		}
	}
}
