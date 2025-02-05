#!/bin/bash
#SBATCH --job-name=hap_align
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni/mummer/align_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/l_dearborni/mummer/align_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30GB
#SBATCH --time=5-00:00:00

#load required modules
module load miniconda3
conda activate mummer

nucmer --coords -p l_dearborni_h1tg000019l_h2tg000034l /hb/home/omoosman/owen/zoarcoidei/data/assemblies/l_dearborni/l_dearborni_h1tg000019l.fasta /hb/home/omoosman/owen/zoarcoidei/data/assemblies/l_dearborni/l_dearborni_h2tg000034l.fasta

