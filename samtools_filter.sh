#!/bin/bash
#SBATCH --job-name=zoarcoidei_filt
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/filt_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/filt_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --time=1-00:00:00


#load module

module load samtools

#filter with samtools and index output

samtools view -bh /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/a_lupus/a_lupus_alignment.sorted.bam h1tg000009l h1tg000038l h1tg000066l h2tg000001l h2tg000051l h2tg000193l -o /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/a_lupus/a_lupus_filt_AFP.bam

samtools index /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/a_lupus/a_lupus_filt_AFP.bam

exonerate --model protein2genome \ 
  --query /path/to/Mamericanus_AFP.txt \
  --target /path/to/pgunn/hap_2_genome.fa \
  --showtargetgff TRUE \
  --showquerygff FALSE \
  --minintron 0 \
  --maxintron 10000 \
  --showalignment FALSE --showcigar FALSE \
  --ryo "Query: %qi Length: %ql Strand: %qs Target: %ti Range: %tcb-%tce\n" \
  > /path/to/output/pgunn_hap2_AFP_annotations.gff
