#!/bin/bash
#SBATCH --job-name=zoarcoidei_exonerate
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/exonerate/query_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/exonerate/query_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --time=0-1:00:00

#load required modules

module load miniconda3

#define variables

config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/exonerate/$name

# Create the output directory 

mkdir -p "$out"

#convert gfa to fasta

conda activate gfatools

gfatools gfa2fa "$in/${name}_hifi.asm.bp.hap1.p_ctg.gfa" > "$in/${name}_hap1_ctg.fasta" 

gfatools gfa2fa "$in/${name}_hifi.asm.bp.hap2.p_ctg.gfa" > "$in/${name}_hap2_ctg.fasta"

gfatools gfa2fa "$in/${name}_hifi.asm.bp.r_utg.gfa" > "$in/${name}_r_utg.fasta"

conda deactivate

#exonerate query for AFPs 

conda activate exonerate

exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/analysis/Mamericanus_AFP.txt --target "$in/${name}_hap1_ctg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/${name}_hap1_AFP_exonerate_output.txt"

exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/analysis/Mamericanus_AFP.txt --target "$in/${name}_hap2_ctg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/${name}_hap2_AFP_exonerate_output.txt"

exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/analysis/Mamericanus_AFP.txt --target "$in/${name}_r_utg.fasta"  --showtargetgff TRUE --showquerygff FALSE > "$out/${name}_r_utg_AFP_exonerate_output.txt"
