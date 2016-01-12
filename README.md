# osmotic

1: run trimmomatic

	15.02.17.trimmomatic_array.sh 

2: run bwa mem, merge and index bams

	bwa_mem_paired_array.sh
	bwa_mem_unpaired_array.sh
	mergebams.sh
	indexbam.sh

3: run featurecounts

	15.11.13.featurecounts.sh

4: analyze counts in edgeR

	osmotic_edgeR_script.R
