#!/bin/bash
#SBATCH --job-name=zoarcoidei_flanking_gene
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/exonerate/flank_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/exonerate/flank_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --time=0-1:00:00

#load required modules

module load miniconda3
conda activate exonerate

#define variables

config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/exonerate/$name

# Create the output directory 

mkdir -p "$out"

#exonerate query for haps 1 and 2
exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/data/Gacul_ldb3b.fa \
  --target "$in/${name}_ref.fasta" \
  --showtargetgff TRUE \
  --showquerygff FALSE \
  --showalignment TRUE --showcigar FALSE \
  --ryo "Query: %qi Length: %ql Strand: %qs Target: %ti Range: %tcb-%tce\n" \
  > "$out/${name}_ldb3b.gff"


  exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/data/Gacul_sncgb.fa \
  --target "$in/${name}_ref.fasta" \
  --showtargetgff TRUE \
  --showquerygff FALSE \
  --showalignment TRUE --showcigar FALSE \
  --ryo "Query: %qi Length: %ql Strand: %qs Target: %ti Range: %tcb-%tce\n" \
  > "$out/${name}_sncgb.gff"


#exonerate quary for raw utg
  exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/data/Gacul_ldb3b.fa \
  --target "$in/${name}_r_utg.fasta" \
  --showtargetgff TRUE \
  --showquerygff FALSE \
  --showalignment TRUE --showcigar FALSE \
  --ryo "Query: %qi Length: %ql Strand: %qs Target: %ti Range: %tcb-%tce\n" \
  > "$out/${name}_ldb3b.gff"


  exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/data/Gacul_sncgb.fa \
  --target "$in/${name}_r_utg.fasta" \
  --showtargetgff TRUE \
  --showquerygff FALSE \
  --showalignment TRUE --showcigar FALSE \
  --ryo "Query: %qi Length: %ql Strand: %qs Target: %ti Range: %tcb-%tce\n" \
  > "$out/${name}_sncgb.gff"
