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

module load samtools
module load miniconda3
conda activate pbmm2 

# defining variables

name=l_dearborni
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name
file=/hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz


#set the tmpdir

export TMPDIR=/hb/scratch/omoosman


#check proportion of reads that aligned in orginal alignment 

samtools view "$out/${name}_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f read_ids.txt aligned_reads.txt | wc -l

## alignment with ccs preset

pbmm2 align --sort --preset CCS "$in/${name}_ref.fasta" $file "$out/${name}_ccs_alignment.sorted.bam"

#filter for only regions with AFPs

samtools view -bh "$out/${name}_ccs_alignment.sorted.bam" h1tg000019l h2tg000034l -o "$out/${name}_ccs_alignment.filt.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_ccs_alignment.filt.sorted.bam" > "$out/${name}_ccs_m_a_alignment.filt.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_ccs_alignment.filt.sorted.bam" > "$out/${name}_ccs_p_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_ccs_m_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_ccs_p_a_alignment.filt.sorted.bam"



## alignment with less stringent parameters

pbmm2 align --sort --preset CCS --min-idt 0.7 --min-score 50 --best-n 5 "$in/${name}_ref.fasta" $file "$out/${name}_lenient_alignment.sorted.bam"

#filter for only regions with AFPs

samtools view -bh "$out/${name}_lenient_alignment.sorted.bam" h1tg000019l h2tg000034l -o "$out/${name}_lenient_alignment.filt.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_lenient_alignment.filt.sorted.bam" > "$out/${name}_lenient_m_a_alignment.filt.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_lenient_alignment.filt.sorted.bam" > "$out/${name}_lenient_p_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_lenient_m_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_lenient_p_a_alignment.filt.sorted.bam"
