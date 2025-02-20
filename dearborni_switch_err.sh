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

name=l_dearborni
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name
file=/hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz


#set the tmpdir

export TMPDIR=/hb/scratch/omoosman


# alignment with ccs preset

pbmm2 align --sort --preset CCS "$in/${name}_ref.fasta" $file "$out/${name}_ccs_alignment.sorted.bam"

# alignment with less stringent parameters

pbmm2 align --sort --preset CCS --min-idt 0.7 --min-score 50 --best-n 5 "$in/${name}_ref.fasta" $file "$out/${name}_lenient_alignment.sorted.bam"
