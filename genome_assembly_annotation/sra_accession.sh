#!/bin/bash
#SBATCH --job-name=SRA_accession
#SBATCH --time=0-24:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --output=/hb/groups/kelley_training/owen/zoarcoidei/data/err_out/sra-%A_%a.out # output file
#SBATCH --error=/hb/groups/kelley_training/owen/zoarcoidei/data/err_out/sra-%A_%a.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=100GB

#define SRA accesion ID & species name
config_file=/hb/groups/kelley_training/owen/zoarcoidei/data/sra_id.txt
id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $1}')
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$config_file" | awk '{print $2}')

#load required modules
module load sratoolkit

#make required directory and move to it
cd /hb/scratch/$USER 

#prefetch accession into scratch
prefetch "$id" --max-size 100g  -O /hb/scratch/$USER 

#convert to fastq from scratch directory 
fasterq-dump "$id"   

#compress and copy files to correct directory
gzip /hb/scratch/$USER/*.fastq
mkdir -p /hb/home/omoosman/owen/zoarcoidei/data/raw_hifi/$name
cp /hb/scratch/$USER/*.fastq.gz /hb/home/omoosman/owen/zoarcoidei/data/raw_hifi/$name
rm *.fastq.gz

