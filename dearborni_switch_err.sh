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
conda activate pbmm2 

# defining variables

name=l_dearborni
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name
file=/hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz


#set the tmpdir

export TMPDIR=/hb/scratch/omoosman

###original alignment 

##check proportion of reads that aligned in orginal alignment 

#calculate number of reads of each type
m=$(wc -l "$out/${name}_m.txt" | awk '{print $1}')
a=$(wc -l "$out/${name}_a.txt" | awk '{print $1}')
p=$(wc -l "$out/${name}_p.txt" | awk '{print $1}')

#calculate the number of reads of each type in the .bam file
o_m=$(samtools view "$out/${name}_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_m.txt" | wc -l)
o_a=$(samtools view "$out/${name}_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_a.txt" | wc -l)
o_p=$(samtools view "$out/${name}_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_p.txt" | wc -l)

#calculate proportions with decimal places
proportion_o_m=$(awk "BEGIN {print $o_m / $m}")
proportion_o_a=$(awk "BEGIN {print $o_a / $a}")
proportion_o_p=$(awk "BEGIN {print $o_p / $p}")

#output in easily readable format
echo "The proportion of hap1 reads that align to the original alignment is $proportion_o_p" >> "$out/proportion.txt"
echo "The proportion of hap2 reads that align to the original alignment is $proportion_o_m" >> "$out/proportion.txt"
echo "The proportion of unphased reads that align to the original alignment is $proportion_o_a" >> "$out/proportion.txt"
echo "" >> "$out/proportion.txt"


###alignment with ccs preset

pbmm2 align --sort --preset CCS "$in/${name}_ref.fasta" $file "$out/${name}_ccs_alignment.sorted.bam"

#filter for only regions with AFPs

samtools view -bh "$out/${name}_ccs_alignment.sorted.bam" h1tg000019l h2tg000034l -o "$out/${name}_ccs_alignment.filt.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_ccs_alignment.filt.sorted.bam" > "$out/${name}_ccs_m_a_alignment.filt.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_ccs_alignment.filt.sorted.bam" > "$out/${name}_ccs_p_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_ccs_m_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_ccs_p_a_alignment.filt.sorted.bam"

##check proportion of reads that aligned in ccs alignment 

#calculate the number of reads of each type in the .bam file
ccs_m=$(samtools view "$out/${name}_ccs_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_m.txt" | wc -l)
ccs_a=$(samtools view "$out/${name}_ccs_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_a.txt" | wc -l)
ccs_p=$(samtools view "$out/${name}_ccs_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_p.txt" | wc -l)

#calculate proportions with decimal places
proportion_ccs_m=$(awk "BEGIN {print $ccs_m / $m}")
proportion_ccs_a=$(awk "BEGIN {print $ccs_a / $a}")
proportion_ccs_p=$(awk "BEGIN {print $ccs_p / $p}")

#output in easily readable format
echo "The proportion of hap1 reads that align to the ccs aligment is $proportion_ccs_p" >> "$out/proportion.txt"
echo "The proportion of hap2 reads that align to the ccs aligment is $proportion_ccs_m" >> "$out/proportion.txt"
echo "The proportion of unphased reads that align to the ccs aligment is $proportion_ccs_a" >> "$out/proportion.txt"
echo "" >> "$out/proportion.txt"


###alignment with less stringent parameters

pbmm2 align --sort --preset CCS -min-id-perc 0.7 -min-concordance-perc 50 --best-n 5 "$in/${name}_ref.fasta" $file "$out/${name}_lenient_alignment.sorted.bam"

#filter for only regions with AFPs

samtools view -bh "$out/${name}_lenient_alignment.sorted.bam" h1tg000019l h2tg000034l -o "$out/${name}_lenient_alignment.filt.sorted.bam"

#filter by reads aligning to each haplotype and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_lenient_alignment.filt.sorted.bam" > "$out/${name}_lenient_m_a_alignment.filt.sorted.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_lenient_alignment.filt.sorted.bam" > "$out/${name}_lenient_p_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_lenient_m_a_alignment.filt.sorted.bam"

samtools index "$out/${name}_lenient_p_a_alignment.filt.sorted.bam"

##check proportion of reads that aligned in lenient alignment 

#calculate the number of reads of each type in the .bam file
lax_m=$(samtools view "$out/${name}_lenient_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_m.txt" | wc -l)
lax_a=$(samtools view "$out/${name}_lenient_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_a.txt" | wc -l)
lax_p=$(samtools view "$out/${name}_lenient_alignment.sorted.bam" | awk '{print $1}' | grep -F -f "$out/${name}_p.txt" | wc -l)

#calculate proportions with decimal places
proportion_lax_m=$(awk "BEGIN {print $lax_m / $m}")
proportion_lax_a=$(awk "BEGIN {print $lax_a / $a}")
proportion_lax_p=$(awk "BEGIN {print $lax_p / $p}")

#output in easily readable format
echo "The proportion of hap1 reads that align to the lenient aligment is $proportion_lax_p" >> "$out/proportion.txt"
echo "The proportion of hap2 reads that align to the lenient aligment is $proportion_lax_m" >> "$out/proportion.txt"
echo "The proportion of unphased reads that align to the lenient aligment is $proportion_lax_a" >> "$out/proportion.txt"
echo "" >> "$out/proportion.txt"
