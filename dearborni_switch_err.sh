#!/bin/bash
#SBATCH --job-name=zoarcoidei_switch_err
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/switch_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/switch_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=70GB
#SBATCH --time=5-00:00:00


###### Realignment of haploptype specific reads to haplotype assemblies for dearborni, using different parameters #######

# loading required modules + conda environment

module load miniconda3
conda activate pbmm2 

# defining variables

in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/l_dearborni
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni
