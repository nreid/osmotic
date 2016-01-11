#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=1-130
#SBATCH -p hi
#SBATCH --cpus-per-task=2


#load software for this job
module load boost/1.55.0
module load samtools/1.2
module load bowtie2/2.2.5
module load tophat/2.1.0
module load bwa/0.7.9a



cd ~/rnaseq/alignments/merged

file=$(ls *bam | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)

samtools index $file
