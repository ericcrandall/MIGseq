#trim the first 14 bases (the msat + anchor sequences of PCR1 primers) from Read 2, and the last 46 bases to leave a fragment of 90bp
#trim 60 bases from the 3' end of read 1 to leave a fragment of 90bp
for i in *R2*.fastq

	do 
		cat $i | fastx_trimmer -Q33 -f 15 -o temp.fastq
		cat temp.fastq | fastx_trimmer -Q33 -t 46 -o ../2trimmed/$(basename $i .fastq)_trim.fastq
	done

for i in *R1*.fastq
	do
		cat $i | fastx_trimmer -Q33 -t 60 -o ../2trimmed/$(basename $i .fastq)_trim.fastq 
	done
