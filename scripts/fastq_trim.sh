for i in *R2*.fastq

	do 
		cat $i | fastx_trimmer -Q33 -t 14 -o ../trimmed/$(basename $i .fastq)_fqf.fastq
	done 
