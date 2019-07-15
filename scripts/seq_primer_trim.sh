#!/bin/bash

for i in *R1*.fastq

	do 
		tagdust -ref /projects/crandall/MIGseq/Read2Seq.txt -dust 100 -1 R:N $i  -o ../4dusted/$(basename $i .fastq)_dust
	done 


for i in *R2*.fastq

	do 
		tagdust -ref /projects/crandall/MIGseq/Read1Seq.txt -dust 100 -1 R:N $i -o ../4dusted/$(basename $i .fastq)_dust
	done 
