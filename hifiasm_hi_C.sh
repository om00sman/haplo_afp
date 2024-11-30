#!/bin/bash
#SBATCH --job-name=zoarcoidei_hifiasm
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/data/err_out/hifiasm_onlyhifi_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/data/err_out/hifiasm_onlyhifi_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --time=5-00:00:00
#SBATCH --partition=256x44

# Load dependencies
module load miniconda3
conda activate hifiasm

#defining species name and path to fastq.gz
config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')
id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $1}')

#define input and output files
input=$(ls /hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz)
h1=/hb/home/omoosman/owen/zoarcoidei/analysis/trimmed_hi_c/$name/"${id}_1_val_1.fq.gz" 
h2=/hb/home/omoosman/owen/zoarcoidei/analysis/trimmed_hi_c/$name/"${id}_2_val_2.fq.gz" 
output=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name

# Make out directory and move to directory with input files
mkdir -p "$output"
cd /hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name

# hifi reads only
hifiasm -o "${output}/${name}_hifi.asm" -t20 --write-paf --h1 "$h1" --h2 "$h2" "$input" 2> "${output}/${name}_hifi.asm.log"
