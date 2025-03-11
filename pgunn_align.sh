#!/bin/bash
#SBATCH --job-name=dearborni_utg_align
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/dearborni_utg_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/dearborni_utg_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50GB
#SBATCH --time=5-00:00:00


###### Realignment of hap1 reads to simulated switch error #####

# loading required modules + conda environment

module load samtools
module load miniconda3
conda activate pbmm2 

# defining variables

name=p_gunnellus
in=/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/$name
out=/hb/groups/kelley_training/owen/zoarcoidei/analysis/realignment/$name


#set the tmpdir

export TMPDIR=/hb/scratch/$USER

#pull hap1 reads from file

seqkit grep -f "$out/${name}_p_a.txt" "$out/${name}_combined.fastq.gz" -o "$out/${name}_hap1_reads.fastq.gz"

###alignment to simulated switch err

pbmm2 align --sort "$out/pgunn_sim_sw_err.fasta" "$out/${name}_hap1_reads.fastq.gz" "$out/${name}_hap1_sim_sw_err_alignment.sorted.bam"

#index output

samtools index "$out/${name}_hap1_sim_sw_err_alignment.sorted.bam"


