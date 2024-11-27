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
#SBATCH --mem=35GB
#SBATCH --time=1-00:00:00


#load module

module load samtools

#filter with samtools

samtools view -bh /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/p_gunnellus/p_gunnellus_alignment.sorted.bam h1tg000157l h2tg000025l -o /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/p_gunnellus/p_gunnellus_filt_AFP.bam
