#!/bin/bash
#SBATCH --job-name=zoarcoidei_filt2
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/filt2_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/filt2_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --time=1-00:00:00


#load module

module load samtools

# defining variables

config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name



#filter with samtools and index output

samtools view -bh /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/a_lupus/a_lupus_alignment.sorted.bam h1tg000009l h1tg000038l h1tg000066l h2tg000001l h2tg000051l h2tg000193l -o /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/a_lupus/a_lupus_filt_AFP.bam

samtools index /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/a_lupus/a_lupus_filt_AFP.bam
