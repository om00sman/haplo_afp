#!/bin/bash
#SBATCH --job-name=hap_align
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/hap_align_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/hap_align_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=70GB
#SBATCH --time=5-00:00:00

#load required modules
module load minimap2
module load samtools

#define variables
config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name

#perform alignment

minimap2 -t 12 -ax asm5 "$in/${name}_hap1_ctg.fasta" "$in/${name}_hap2_ctg.fasta" > "$out/${name}_haptohap.sam"

#convert output to BAM and sort and index
samtools view -b "$out/${name}_haptohap.sam" > "$out/${name}_haptohap.bam"
samtools sort "$out/${name}_haptohap.bam" -o "$out/${name}_haptohap.sorted.bam"
samtools index "$out/${name}_haptohap.sorted.bam"
