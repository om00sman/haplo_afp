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

#define variables

config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=z_americanus
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name

#filter with samtools and index output
"$out/${name}_filt_AFP_p_a.bam"

samtools view -bh /hb/home/omoosman/owen/zoarcoidei/analysis/realignment/z_americanus/z_americanus_alignment.sorted.bam h1tg000092l h1tg000766l h2tg000180l h2tg000104l -o "$out/${name}_filt_AFP_p_a.bam"

samtools index "$out/${name}_filt_AFP_p_a.bam"
