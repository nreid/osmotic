#!/bin/bash
#SBATCH -p hi
#SBATCH --cpus-per-task=24
#SBATCH --mem=16000
###### number of nodes
###### number of processors

module load samtools/1.2 
module load boost/1.55.0

#make sure temp files are written to scratch by running script from there. 
cd /scratch/
mkdir nreid
cd nreid

#various relevant directories and files
genomebase=/home/nreid/popgen/kfish3/killifish20130322asm
genome=/home/nreid/popgen/kfish3/killifish20130322asm.fa

fc=/home/nreid/bin/subread-1.4.6-p4-source/bin/featureCounts
saf=/home/nreid/popgen/kfish3/kfish2rae5g.main.pub.corrected.exons.saf
inbam=$(cat /home/nreid/rnaseq/bwa.bams.list | tr '\n' ',' | sed 's/,/ /g')

root=osmotic.counts
outdir=/home/nreid/rnaseq/featurecounts/$root.fc
outdir2=/home/nreid/rnaseq/featurecounts/$root.meta.fc

$fc -Q 10 -T 24    -p -a $saf -F SAF -o $outdir2 $inbam
$fc -Q 10 -T 24 -f    -a $saf -F SAF -o $outdir  $inbam

cd ..
rm -r nreid
