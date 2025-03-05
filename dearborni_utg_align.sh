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

###alignment to hap1 shortest utg

pbmm2 align --sort "$in/ldear__hap1_shortest.fasta" $file "$out/${name}_hap1s_alignment.sorted.bam"

#index output

samtools index "$out/${name}_hap1s_alignment.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_hap1s_alignment.sorted.bam" > "$out/${name}_hap1s_m_a_alignment.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_hap1s_alignment.sorted.bam" > "$out/${name}_hap1s_p_a_alignment.sorted.bam"

samtools index "$out/${name}_hap1s_m_a_alignment.sorted.bam"

samtools index "$out/${name}_hap1s_p_a_alignment.sorted.bam"

##check proportion of reads that aligned in utg alignment 

#calculate number of reads of each type
m=$(wc -l "$out/${name}_m.txt" | awk '{print $1}')
a=$(wc -l "$out/${name}_a.txt" | awk '{print $1}')
p=$(wc -l "$out/${name}_p.txt" | awk '{print $1}')

#calculate the number of reads of each type in the .bam file
hap1s_m=$(samtools view "$out/${name}_hap1s_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_m.txt" | wc -l)
hap1s_a=$(samtools view "$out/${name}_hap1s_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_a.txt" | wc -l)
hap1s_p=$(samtools view "$out/${name}_hap1s_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_p.txt" | wc -l)

#calculate proportions with decimal places
proportion_1s_m=$(awk "BEGIN {print $hap1s_m / $m}")
proportion_1s_a=$(awk "BEGIN {print $hap1s_a / $a}")
proportion_1s_p=$(awk "BEGIN {print $hap1s_p / $p}")

#output in easily readable format
echo "The proportion of hap1 reads that align to the hap1 shortest utg assembly is $proportion_1s_p" >> "$out/proportion.txt"
echo "The proportion of hap2 reads that align to the hap1 shortest utg assembly is $proportion_1s_m" >> "$out/proportion.txt"
echo "The proportion of unphased reads that align to the hap1 shortest utg assembly is $proportion_1s_a" >> "$out/proportion.txt"
echo "" >> "$out/proportion.txt"

###alignment to hap1 longest utg

pbmm2 align --sort "$in/ldear__hap1_longest.fasta" $file "$out/${name}_hap1l_alignment.sorted.bam"

#index output

samtools index "$out/${name}_hap1l_alignment.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_hap1l_alignment.sorted.bam" > "$out/${name}_hap1l_m_a_alignment.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_hap1l_alignment.sorted.bam" > "$out/${name}_hap1l_p_a_alignment.sorted.bam"

samtools index "$out/${name}_hap1l_m_a_alignment.sorted.bam"

samtools index "$out/${name}_hap1l_p_a_alignment.sorted.bam"

##check proportion of reads that aligned in utg alignment 

#calculate number of reads of each type
m=$(wc -l "$out/${name}_m.txt" | awk '{print $1}')
a=$(wc -l "$out/${name}_a.txt" | awk '{print $1}')
p=$(wc -l "$out/${name}_p.txt" | awk '{print $1}')

#calculate the number of reads of each type in the .bam file
hap1l_m=$(samtools view "$out/${name}_hap1l_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_m.txt" | wc -l)
hap1l_a=$(samtools view "$out/${name}_hap1l_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_a.txt" | wc -l)
hap1l_p=$(samtools view "$out/${name}_hap1l_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_p.txt" | wc -l)

#calculate proportions with decimal places
proportion_1l_m=$(awk "BEGIN {print $hap1l_m / $m}")
proportion_1l_a=$(awk "BEGIN {print $hap1l_a / $a}")
proportion_1l_p=$(awk "BEGIN {print $hap1l_p / $p}")

#output in easily readable format
echo "The proportion of hap1 reads that align to the hap1 longest utg assembly is $proportion_1l_p" >> "$out/proportion.txt"
echo "The proportion of hap2 reads that align to the hap1 longest utg assembly is $proportion_1l_m" >> "$out/proportion.txt"
echo "The proportion of unphased reads that align to the hap1 longest utg assembly is $proportion_1l_a" >> "$out/proportion.txt"
echo "" >> "$out/proportion.txt"

###alignment to hap2 shortest utg

pbmm2 align --sort "$in/ldear__hap2_shortest.fasta" $file "$out/${name}_hap2s_alignment.sorted.bam"

#index output

samtools index "$out/${name}_hap2s_alignment.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_hap2s_alignment.sorted.bam" > "$out/${name}_hap2s_m_a_alignment.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_hap2s_alignment.sorted.bam" > "$out/${name}_hap2s_p_a_alignment.sorted.bam"

samtools index "$out/${name}_hap2s_m_a_alignment.sorted.bam"

samtools index "$out/${name}_hap2s_p_a_alignment.sorted.bam"

##check proportion of reads that aligned in utg alignment 

#calculate number of reads of each type
m=$(wc -l "$out/${name}_m.txt" | awk '{print $1}')
a=$(wc -l "$out/${name}_a.txt" | awk '{print $1}')
p=$(wc -l "$out/${name}_p.txt" | awk '{print $1}')

#calculate the number of reads of each type in the .bam file
hap2s_m=$(samtools view "$out/${name}_hap2s_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_m.txt" | wc -l)
hap2s_a=$(samtools view "$out/${name}_hap2s_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_a.txt" | wc -l)
hap2s_p=$(samtools view "$out/${name}_hap2s_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_p.txt" | wc -l)

#calculate proportions with decimal places
proportion_2s_m=$(awk "BEGIN {print $hap2s_m / $m}")
proportion_2s_a=$(awk "BEGIN {print $hap2s_a / $a}")
proportion_2s_p=$(awk "BEGIN {print $hap2s_p / $p}")

#output in easily readable format
echo "The proportion of hap1 reads that align to the hap2 shortest utg assembly is $proportion_2s_p" >> "$out/proportion.txt"
echo "The proportion of hap2 reads that align to the hap2 shortest utg assembly is $proportion_2s_m" >> "$out/proportion.txt"
echo "The proportion of unphased reads that align to the hap2 shortest utg assembly is $proportion_2s_a" >> "$out/proportion.txt"
echo "" >> "$out/proportion.txt"

###alignment to hap2 longest utg

pbmm2 align --sort "$in/ldear__hap2_longest.fasta" $file "$out/${name}_hap2l_alignment.sorted.bam"

#index output

samtools index "$out/${name}_hap2l_alignment.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_hap2l_alignment.sorted.bam" > "$out/${name}_hap2l_m_a_alignment.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_hap2l_alignment.sorted.bam" > "$out/${name}_hap2l_p_a_alignment.sorted.bam"

samtools index "$out/${name}_hap2l_m_a_alignment.sorted.bam"

samtools index "$out/${name}_hap2l_p_a_alignment.sorted.bam"

##check proportion of reads that aligned in utg alignment 

#calculate number of reads of each type
m=$(wc -l "$out/${name}_m.txt" | awk '{print $1}')
a=$(wc -l "$out/${name}_a.txt" | awk '{print $1}')
p=$(wc -l "$out/${name}_p.txt" | awk '{print $1}')

#calculate the number of reads of each type in the .bam file
hap2l_m=$(samtools view "$out/${name}_hap2l_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_m.txt" | wc -l)
hap2l_a=$(samtools view "$out/${name}_hap2l_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_a.txt" | wc -l)
hap2l_p=$(samtools view "$out/${name}_hap2l_alignment.sorted.bam" | awk '{print $1}' | sort | uniq | grep -F -f "$out/${name}_p.txt" | wc -l)

#calculate proportions with decimal places
proportion_2l_m=$(awk "BEGIN {print $hap2l_m / $m}")
proportion_2l_a=$(awk "BEGIN {print $hap2l_a / $a}")
proportion_2l_p=$(awk "BEGIN {print $hap2l_p / $p}")

#output in easily readable format
echo "The proportion of hap1 reads that align to the hap2 longest utg assembly is $proportion_2l_p" >> "$out/proportion.txt"
echo "The proportion of hap2 reads that align to the hap2 longest utg assembly is $proportion_2l_m" >> "$out/proportion.txt"
echo "The proportion of unphased reads that align to the hap2 longest utg assembly is $proportion_2l_a" >> "$out/proportion.txt"
echo "" >> "$out/proportion.txt"

