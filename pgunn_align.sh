#!/bin/bash
#SBATCH --job-name=dearborni_utg_align
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/dearborni_utg_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/dearborni_utg_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50GB
#SBATCH --time=5-00:00:00


###### Realignment of hap1 reads to simulated switch error #####

# loading required modules + conda environment

module load samtools
module load miniconda3
conda activate pbmm2 

# defining variables

name=p_gunnellus
in=/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/$name
out=/hb/groups/kelley_training/owen/zoarcoidei/analysis/realignment/$name


#set the tmpdir

export TMPDIR=/hb/scratch/$USER

#pull hap1 reads from file

seqkit grep -f  seqs.fq.gz -o result.fq.gz


###alignment to hap1 shortest utg

pbmm2 align --sort "$in/ldear__hap1_shortest.fasta" $file "$out/${name}_hap1s_alignment.sorted.bam"

#index output

samtools index "$out/${name}_hap1s_alignment.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_hap1s_alignment.sorted.bam" > "$out/${name}_hap1s_m_a_alignment.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_hap1s_alignment.sorted.bam" > "$out/${name}_hap1s_p_a_alignment.sorted.bam"

samtools index "$out/${name}_hap1s_m_a_alignment.sorted.bam"

samtools index "$out/${name}_hap1s_p_a_alignment.sorted.bam"

