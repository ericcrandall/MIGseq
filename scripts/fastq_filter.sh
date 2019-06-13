for i in *.fastq

	do 
		cat $i | fastq_quality_filter -Q33 -q 30 -p 20 -o ../filtered/$(basename $i .fastq)_fqf.fastq
	done 
