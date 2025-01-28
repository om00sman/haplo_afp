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

#define variables
minimap2 -t 12 -ax asm5 ref.fa asm.fa > aln.sam
