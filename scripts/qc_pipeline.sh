#!/bin/bash

#create directory structure
mkdir 1original
mkdir 2trimmed
mkdir 3filtered
mkdir 4dusted
mkdir 5concatenated
mkdir 6preprocessed
mkdir 7stacks

#move everything into original
mv *.fastq ./1original

#Trim
echo "Running fastx_trimmer to remove 14 msat bases from read 2 and truncate all reads to 90bp"
cd 1original
../../fastq_trim.sh
cd ..

#Filter
echo "Running fastq_quality_filter to trim low quality bases and remove low quality reads"
cd 2trimmed
../../fastq_filter.sh
cd ..

#Tagdust - remove sequencing primers and poly A sequences
echo "Running Tagdust to remove sequencing primers and low complexity reads"
cd 3filtered
../../seq_primer_trim.sh
cd ..

# Concatenate the read 1 and read 2 together for each individual (as suggested by Suyama et al. 2015)
echo "Concatenate read 1 and read 2 and rename to just the sample name"
#rename to just the individual identifier for input into stacks
cd 4dusted
#move the extraneous files out of the way
mkdir tagdust_logs
mv *logfile* ./tagdust_logs
mkdir unextracted
mv *un.fq ./unextracted
#this for loop finds all read 1s, takes the sample name (delimited by -) and
#concatenates the two files with this prefix, placing the resulting file, named by the sample name, in the final directory
for filename in *R1*dust.fq 
    do
      prefix=$(echo $filename | sed 's/^\([a-z]*-[0-9]*-[0-9]*\).*/\1/')
      echo "now concatenating $prefix"
      cat $prefix* > ../5concatenated/$prefix.fq
    done
cd ..

# remove sequences with Ns and zip them up
echo "Running process_shortreads to remove reads with N and zip final files"
process_shortreads -p ./5concatenated -o ./6preprocessed -c

cd ..
