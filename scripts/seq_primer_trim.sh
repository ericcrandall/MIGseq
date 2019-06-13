#!/bin/bash

for i in *R1*.fastq

        do
          	tagdust -ref /projects/crandall/MIGseq/Read2Seq.txt -dust 100 -1 R:N $i  -o ../seq_primer_trimmed/$(basename $i .fastq)_final.fastq
        done


for i in *R2*.fastq

        do
          	tagdust -ref /projects/crandall/MIGseq/Read1Seq.txt -dust 100 -1 R:N $i -o ../seq_primer_trimmed/$(basename $i .fastq)_final.fastq
        done