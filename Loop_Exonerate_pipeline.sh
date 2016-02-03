#!/usr/bin/bash

########################################################################################
#shell script to run post aTRAM pipeline created by J M Allen.
#script written by B M Boyd, UIUC
#Date created: 6-May-2014
_______________________________________________________________________________________
#USAGE:
#Usage: copy this shell script into your directory with best files 
#Usage: use emacs to input variables in block below
#Usage: sh FIXED_copy_and_run_whole_pipeline_at_once.sh
_______________________________________________________________________________________
#VARIABLES TO DEFINE BEFORE RUNNING (substitute "<<repalce me>>" with your "variable"):
#may replace with "prompt -p var1 var2 narN;" at some point and use echo to list what to input;
#mylib='Stjap';    ### input library name, eg. "Covei"  ### THIS PIPELINE WAS MEANT TO RUN ON A BUNCH OF LIBRARIES IN ONE FOLDER. SO THIS IS CONFUSING
myaTRAMfile='trinity.out.best.ed'; #    I4.abyss.out.best.ed';  ### input rest of file name for aTRAM contigs, eg. "I4.trinity.ed"  ## THIS IS A BIT CONFUSING AS THE FILENAME DOES NOT EXIST UNTIL AFTER THE FIRST SCRIPT and should not include the .fasta.
summaryfile='SummaryIschnocera8_4_15';  ### input a file name of your choice, this is just a log file 
mypath='/data/ischnocera3/1107Exon_Stitching'; ### pathway to directory, leave off last "/"
myoverlap='10';  ### overlap to use, typically is "10" but can change if needed  #### THIS SHOULD BE 10 FOR THE FIRST ITERATION AND 0 FOR THE SECOND
########################################################################################

#step 1
cp /data/scripts/Exonerate_pipeline/editbestfiles.pl .
echo "step 1, running editbestfiles";
perl editbestfiles.pl;
echo "done with step 1";

#step 2
echo "step 2, running exonerate the first time";
cd /data/target_genes/GET1107ORTH_PEPTIDES/;
for myref in *.pep.fasta; do
######	debugging 5.18.15	
	for lib in  Covei; do    #SpspTapor; do #   Crimm Dobre Famar Fulon Hldiv Ibbis Licap SpspTabor; do
 		mygene=`expr match "$myref" '.*\(PHUM[0-9]*\)'`; 
		 exonerate --model protein2genome /data/target_genes/GET1107ORTH_PEPTIDES/$myref $mypath/$mygene.$lib.$myaTRAMfile.fasta --showvulgar no --showalignment no --ryo "$mygene,$lib,%ql,%qal,%qab,%qae,abyss,%ti\n" >> $mypath/$mygene.results.out;
	done;
	###### need each lib to be done separately and printed to one file then run the fix.csv script.
	 perl /data/scripts/Exonerate_pipeline/fix_csv.pl $mypath/$mygene.results.out $mypath/$mygene.results.csv;
	LC_ALL=C sort -t, -k 1,1d -k 6,6d -k 4,4n $mypath/$mygene.results.csv > $mypath/$mygene.results.sorted.csv;
done;

echo "done with step 2";

#cd $mypath;

#step 3
#cp /data/scripts/Exonerate_pipeline/getcontigs.first.pl .
#echo "step 3, running first get contigs";
#perl getcontigs.first.pl;
#echo "done with step 3";

#step 4
#cp /data/scripts/Exonerate_pipeline/stitch.aTRAM.contigs.BMBversion.pl .
#echo "step 4, stichting";
#perl stitch.aTRAM.contigs.BMBversion.pl $myoverlap $myaTRAMfile $mypath;
#echo "done with step 4";


#step 5
#cp /data/scripts/Exonerate_pipeline/summarystats.first.pl .
#echo "step 5, getting summary data";
#perl summarystats.first.pl >> $summaryfile;
#echo "done with step 5";
#
#step 6
#echo "step 6, running exonerate for second time";
#cd /data/target_genes/GET1107ORTH_PEPTIDES/;
#for secref in *.pep.fasta; do
#  secgene=`expr match "$secref" '.*\(PHUM[0-9]*\)'`;
#  exonerate --model protein2genome /data/target_genes/GET1107ORTH_PEPTIDES/$secref $mypath/$secgene.aTRAM.Exonerate.Round1.OVERLAP.10.fasta  --showvulgar no --showalignment no --ryo "$secgene,%ql,%qal,%qab,%qae,abyss,%ti\n" >> $mypath/$secgene.exonerate2.out;
#  perl /data/scripts/Exonerate_pipeline/fix_csv.pl $mypath/$secgene.exonerate2.out $mypath/$secgene.exonerate2.csv; 
#  LC_ALL=C sort -t, -k 1,1d -k 6,6d -k 4,4n $mypath/$secgene.exonerate2.csv > $mypath/$secgene.exonerate2.csv.sorted;
#done;
#for thiref in *.pep.fasta; do
# thigene=`expr match "$thiref" '.*\(PHUM[0-9]*\)'`;
# exonerate --model protein2genome /data/target_genes/GET1107ORTH_PEPTIDES/$thiref $mypath/$thigene.aTRAM.Exonerate.Round1.OVERLAP.10.fasta  --showvulgar no --showalignment no --verbose 0 --ryo ">$thigene,abyss,%ti,%qab,%qae\n%tcs\n" >> $mypath/$thigene.exonerate2.fasta;
#done;
#cd $mypath;
#cp  /data/scripts/Exonerate_pipeline/editcontigs.pl .
#perl editcontigs.pl;
#echo "done with step 6";

#step 7
#cp /data/scripts/Exonerate_pipeline/getcontigs.second.pl .
#echo "step 7, running second get contigs";
#perl getcontigs.second.pl;
#echo "done with step 7";

#step 8
#cp /data/scripts/Exonerate_pipeline/stitch.Exonerate.contigs.pl .
#echo "step 8, stitching exonerate contigs";
#perl stitch.Exonerate.contigs.pl;
#echo "done with step 8";

#step 9
#cp /data/scripts/Exonerate_pipeline/summarystats.second.pl .
#echo "step 9, generating summary data for the final time";
#perl summarystats.second.pl >> $summaryfile;
#echo "done with step 9";
