#!/bin/bash
#SBATCH --job-name=kalign
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/kalign_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/kalign_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --time=5-00:00:00

#load required modules
module load miniconda3
conda activate kalign

kalign -i /hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/l_dearborni/l_dearborni_contigs.fasta --type dna -o /hb/groups/kelley_training/owen/zoarcoidei/analysis/realignment/l_dearborni/l_dearborni_alignment.fasta

