mkdir filter_logs

for i in *.fastq

	do 
		cat $i | fastq_quality_filter -Q33 -q 30 -p 40 -o ../3filtered/$(basename $i .fastq)_filter.fastq > ./filter_logs/$(basename $i .fastq).txt
	done 
