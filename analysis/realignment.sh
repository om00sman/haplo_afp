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

###### Realignment of haploptype specific reads to haplotype assemblies #######

# loading required modules + conda environment

module load miniconda3
conda activate pbmm2 

# defining variables

config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')
in=/hb/home/omoosman/owen/zoarcoidei/data/assemblies/$name
out=/hb/home/omoosman/owen/zoarcoidei/analysis/realignment/$name
file1=$(ls /hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz | sed -n '1p')
file2=$(ls /hb/groups/kelley_training/owen/zoarcoidei/data/raw_hifi/$name/*.fastq.gz | sed -n '2p')

# Create the output directory 
mkdir -p "$out"

#set the tmpdir

export TMPDIR=/hb/scratch/omoosman/switch_err

# extract read name and phasing information from haplotype specific .gfa files

awk '$9 == "HG:A:m" {print $5}' "$in/${name}_hifi.asm.bp.hap1.p_ctg.noseq.gfa" > "$out/${name}_hap1_m.txt"
awk '$9 == "HG:A:m" {print $5}' "$in/${name}_hifi.asm.bp.hap2.p_ctg.noseq.gfa" > "$out/${name}_hap2_m.txt"
cat "$out/${name}_hap1_m.txt" "$out/${name}_hap2_m.txt" > "$out/${name}_m.txt"

awk '$9 == "HG:A:a" {print $5}' "$in/${name}_hifi.asm.bp.hap1.p_ctg.noseq.gfa" > "$out/${name}_hap1_a.txt"
awk '$9 == "HG:A:a" {print $5}' "$in/${name}_hifi.asm.bp.hap2.p_ctg.noseq.gfa" > "$out/${name}_hap2_a.txt"
cat "$out/${name}_hap1_a.txt" "$out/${name}_hap2_a.txt" > "$out/${name}_a.txt"

awk '$9 == "HG:A:p" {print $5}' "$in/${name}_hifi.asm.bp.hap1.p_ctg.noseq.gfa" > "$out/${name}_hap1_p.txt"
awk '$9 == "HG:A:p" {print $5}' "$in/${name}_hifi.asm.bp.hap2.p_ctg.noseq.gfa" > "$out/${name}_hap2_p.txt"
cat "$out/${name}_hap1_p.txt" "$out/${name}_hap2_p.txt" > "$out/${name}_p.txt"

# combine read haplotype specific read information

cat "$out/${name}_m.txt" "$out/${name}_a.txt" > "$out/${name}_m_a.txt"
cat "$out/${name}_p.txt" "$out/${name}_a.txt" > "$out/${name}_p_a.txt"

#combine and index reference fasta files

cat "$in/${name}_hap1_ctg.fasta" "$in/${name}_hap2_ctg.fasta" > "$in/${name}_ref.fasta"

pbmm2 index "$in/${name}_ref.fasta" "$in/${name}_ref.mmi" 

## alignment for species with two input files

if [[ "$file2" == *.fastq.gz ]]; then

    zcat $file1 $file2 | gzip > "$out/${name}_combined.fastq.gz"

    pbmm2 align "$in/${name}_ref.fasta" "$out/${name}_combined.fastq.gz" "$out/${name}_alignment.sorted.bam" --sort

### alignment for species with one input file

else

    pbmm2 align "$in/${name}_ref.fasta" $file1 "$out/${name}_alignment.sorted.bam" --sort
    
fi


#filter with samtools and index output

samtools view -bh -N "$out/${name}_m_a.txt" "$out/${name}_filt_AFP.bam" > "$out/${name}_filt_AFP_m_a.bam"

samtools view -bh -N "$out/${name}_p_a.txt" "$out/${name}_filt_AFP.bam" > "$out/${name}_filt_AFP_p_a.bam"

samtools index "$out/${name}_filt_AFP_m_a.bam"

samtools index "$out/${name}_filt_AFP_p_a.bam"


#exonerate query of combined fasta

conda activate exonerate

exonerate --model protein2genome --query /hb/home/omoosman/owen/zoarcoidei/analysis/Mamericanus_AFP.txt \
  --target "$in/${name}_ref.fasta" \
  --showtargetgff TRUE \
  --showquerygff FALSE \
  --minintron 0 \
  --maxintron 10000 \
  --showalignment FALSE --showcigar FALSE \
  --ryo "Query: %qi Length: %ql Strand: %qs Target: %ti Range: %tcb-%tce\n" \
  > "$out/${name}_AFP_annotations.gff"
