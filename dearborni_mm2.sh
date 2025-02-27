#!/bin/bash
#SBATCH --job-name=dearborni_switch_err
#SBATCH --mail-user=omoosman@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --output=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/dearborni_%A_%a.out
#SBATCH --error=/hb/home/omoosman/owen/zoarcoidei/analysis/err_out/switch_err/dearborni_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50GB
#SBATCH --time=5-00:00:00


###### Realignment of haploptype specific reads to haplotype assemblies for dearborni, using different parameters #######

# loading required modules + conda environment

module load samtools
module load miniconda3
conda activate mm2

# defining variables

name=l_dearborni
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name
file=/hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz


#set the tmpdir

export TMPDIR=/hb/scratch/omoosman

###alignment with less stringent parameters

minimap2 -ax map-hifi -o "$out/${name}_mm2_alignment.sam" "$in/${name}_ref.fasta" $file

#filter for only regions with AFPs and convert to bam

samtools view -bh "$out/${name}_mm2_alignment.sam" -o "$out/${name}_mm2_alignment.bam"

samtools sort -o "$out/${name}_mm2_alignment.sorted.bam" "$out/${name}_mm2_alignment.bam"

samtools index "$out/${name}_mm2_alignment.sorted.bam"

samtools view -bh "$out/${name}_mm2_alignment.sorted.bam" h1tg000019l h2tg000034l -o "$out/${name}_mm2_alignment.filt.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_mm2_alignment.filt.sorted.bam" > "$out/${name}_mm2_m_a_alignment.filt.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_mm2_alignment.filt.sorted.bam" > "$out/${name}_mm2_p_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_mm2_m_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_mm2_p_a_alignment.filt.sorted.bam"

##check proportion of reads that aligned in mm2 alignment 

#calculate number of reads of each type
m=$(wc -l "$out/${name}_m.txt" | awk '{print $1}')
a=$(wc -l "$out/${name}_a.txt" | awk '{print $1}')
p=$(wc -l "$out/${name}_p.txt" | awk '{print $1}')

#calculate the number of reads of each type in the .bam file
mm2_m=$(samtools view "$out/${name}_mm2_alignment.sorted.bam" | awk '{print $1}' | uniq | grep -F -f "$out/${name}_m.txt" | wc -l)
mm2_a=$(samtools view "$out/${name}_mm2_alignment.sorted.bam" | awk '{print $1}' | uniq | grep -F -f "$out/${name}_a.txt" | wc -l)
mm2_p=$(samtools view "$out/${name}_mm2_alignment.sorted.bam" | awk '{print $1}' | uniq |grep -F -f "$out/${name}_p.txt" | wc -l)

#calculate proportions with decimal places
proportion_mm2_m=$(awk "BEGIN {print $mm2_m / $m}")
proportion_mm2_a=$(awk "BEGIN {print $mm2_a / $a}")
proportion_mm2_p=$(awk "BEGIN {print $mm2_p / $p}")

#output in easily readable format
echo "The proportion of hap1 reads that align to the mm2 aligment is $proportion_mm2_p" >> "$out/proportion.txt"
echo "The proportion of hap2 reads that align to the mm2 aligment is $proportion_mm2_m" >> "$out/proportion.txt"
echo "The proportion of unphased reads that align to the mm2 aligment is $proportion_mm2_a" >> "$out/proportion.txt"
echo "" >> "$out/proportion.txt"
