#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/mafft_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/mafft_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --time=5-00:00:00

#load required modules
module load mafft

"/hb/software/apps/mafft/gnu-7.520/mafftdir/bin/mafft"  --auto --inputorder "l_dearborni_contigs.fasta" > "l_dearborni_contig_aligned.fasta"

