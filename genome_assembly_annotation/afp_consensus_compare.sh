#!/bin/bash
#SBATCH --job-name=consensus_query
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/exonerate/consensus_query_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/exonerate/consensus_query_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=0-1:00:00

#load required modules

module load miniconda3
conda activate exonerate

#define variables

config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/species.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file")
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/exonerate/$name

# Create the output directory 

mkdir -p "$out"

#exonerate query for with Mamericanus_AFP

exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/analysis/Mamericanus_AFP.txt --target "$in/${name}_hap1_ctg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/consensus_compare/${name}_hap1_Mamer_afp.gff"

exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/analysis/Mamericanus_AFP.txt --target "$in/${name}_hap2_ctg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/consensus_compare/${name}_hap2_Mamer_afp.gff"

#exonerate query for consensus AFP

exonerate --model protein2genome --query "$out/${name}_consensus_afp.fasta" --target "$in/${name}_hap1_ctg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/consensus_compare/${name}_hap1_consensus_afp.gff"

exonerate --model protein2genome --query "$out/${name}_consensus_afp.fasta" --target "$in/${name}_hap2_ctg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/consensus_compare/${name}_hap2_consensus_afp.gff"

#exonerate query for raw utg
exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/analysis/Mamericanus_AFP.txt --target "$in/${name}_r_utg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/consensus_compare/${name}_utg_Mamer_afp.gff"

exonerate --model protein2genome --query "$out/${name}_consensus_afp.fasta" --target "$in/${name}_r_utg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/consensus_compare/${name}_utg_consensus_afp.gff"

