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

samtools view -bh /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_maculatus/l_maculatus_alignment.sorted.bam  h1tg000093l h1tg000447c h1tg000021l h1tg000533l h2tg000350l h2tg000090l h2tg000369l -o /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_maculatus/l_maculatus_filt_AFP.bam

samtools index l_maculatus_filt_AFP.bam
