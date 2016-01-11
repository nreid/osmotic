#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=1-1456
#SBATCH -p hi
#SBATCH --cpus-per-task=6
#SBATCH --mem=10000

#load software for this job
module load boost/1.55.0
module load samtools/1.2
module load bowtie2/2.2.5
module load tophat/2.1.0
module load bwa/0.7.9a

#tell bowtie where the indexed genome is
export BOWTIE2_INDEXES=/home/nreid/popgen/kfish3

#various relevant directories and files
genomebase=/home/nreid/popgen/kfish3/killifish20130322asm
genome=/home/nreid/popgen/kfish3/killifish20130322asm.fa
outdir1=/home/nreid/rnaseq/alignments

#unpaired reads 
fq1=$(find ~/rnaseq/ -type f -name "*.gz" | grep ".*R..*ut" | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)

fr=$(echo $fq1 | grep -oP '(?<=L00._)..')

run=$(echo $fq1 | grep -oP '(?<=rawdata/)[^/]+')

lib=$(echo $fq1 | grep -oP '(?<=Sample_)[^_]+')

samp=$(echo $fq1 | grep -oP "(?<=/AWJRDD00._)[^/]+(?=_[ACGT]+)")

lane=$(echo $fq1 | grep -oP "L00.")

sub=$(echo $fq1 | grep -oP "(?<=L00._R._)...")

rg=$(echo \@RG\\tID:$samp.$lib.$run.$lane.$sub\\tPL:Illumina\\tPU:x\\tLB:$lib\\tSM:$samp)

outdir2=$outdir1/$samp

outfile=$samp.$lib.$run.$lane.$fr.$sub.ut.bam

mkdir -p $outdir2

echo $outdir2
echo $outfile

bwa mem -R $rg $genome $fq1 | \
samtools view -b - | \
samtools fixmate -O bam - - | \
samtools sort -O bam -T /scratch/$outfile.temp >$outdir2/$outfile








