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
module load miniconda3
conda activate exonerate

#filter with samtools and index output

samtools view -bh /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni/l_dearborni_alignment.sorted.bam h1tg000072l h1tg000108l h1tg000110l h2tg000007l h2tg000008l h2tg000061l -o /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni/l_dearborni_filt_AFP.bam

samtools index /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni/l_dearborni_filt_AFP.bam

exonerate --model protein2genome \ 
  --query /hb/home/omoosman/owen/zoarcoidei/analysis/Mamericanus_AFP.txt \
  --target /hb/home/omoosman/owen/zoarcoidei/data/assemblies/l_dearborni/l_dearborni_ref.fasta \
  --showtargetgff TRUE \
  --showquerygff FALSE \
  --minintron 0 \
  --maxintron 10000 \
  --showalignment TRUE --showcigar FALSE \
  --ryo "Query: %qi Length: %ql Strand: %qs Target: %ti Range: %tcb-%tce\n" \
  > /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni/l_dearborni_AFP_annotations.gff
