#!/bin/bash
#SBATCH --job-name=zoarcoidei_hifiasm
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/data/err_out/hifiasm_onlyhifi_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/data/err_out/hifiasm_onlyhifi_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80GB
#SBATCH --time=5-00:00:00

# Load dependencies
module load miniconda3
conda activate hifiasm


#defining species name and path to fastq.gz
config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')
out=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
file1=$(ls /hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz | sed -n '1p')
file2=$(ls /hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz | sed -n '2p')

export TMPDIR=/hb/scratch/omoosman/switch_err

# Make out directory and move to directory with input files
mkdir -p "$out"

## for species with two input files

if [[ "$file2" == *.fastq.gz ]]; then

  #assembly
  hifiasm -o "${out}/${name}_hifi.asm" -t20 --write-paf "$file1" "$file2" 2> "${out}/${name}_hifi.asm.log"

### for species with one input file

else

  #assembly
  hifiasm -o "${out}/${name}_hifi.asm" -t20 --write-paf "$file1" 2> "${out}/${name}_hifi.asm.log"
fi
