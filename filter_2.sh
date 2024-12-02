#!/bin/bash
#SBATCH --job-name=zoarcoidei_filt2
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/filt2_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/filt2_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --time=1-00:00:00


#load module

module load samtools

# defining variables

config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name

# combine read haplotype specific read information

cat "$out/${name}_m.txt" "$out/${name}_a.txt" > "$out/${name}_m_a.txt"
cat "$out/${name}_p.txt" "$out/${name}_a.txt" > "$out/${name}_p_a.txt"

#filter with samtools and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_filt_AFP.bam" > "$out/${name}_filt_AFP_m_a.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_filt_AFP.bam" > "$out/${name}_filt_AFP_p_a.bam"

samtools index "$out/${name}_filt_AFP_m_a.bam"

samtools index "$out/${name}_filt_AFP_p_a.bam"
