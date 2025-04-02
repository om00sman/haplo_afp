# haplo_afp

## Authors

Owen Moosman, University of California, Santa Cruz, Department of Molecular, Cell, & Developmental Biology

Samuel Bogan, University of California, Santa Cruz, Department of Ecology and Evolutionary Biology

Joanna Kelley, University of California, Santa Cruz, Department of Ecology and Evolutionary Biology

## Description

This repository contains code used by Owen Moosman, Samuel Bogan, and Joanna Kelley for analyzing copy number variation of antifreeze protein (AFP) genes between haplotypes of polar fish in the suborder Zoarcoidei. Bioinformatic tools developed and used for this analysis can be found elsewhere. *gfa_parser* is hosted and described at the Github repository https://github.com/snbogan/gfa_parser. *switch_error_screen* is hosted and described at the Github repository https://github.com/snbogan/switch_error_screen. All R scripts were run using R version 4.3.2. 

This research was funded by the United States National Science Foundation, Office of Polar Programs (Award number 2312253) 

## Data Availabilty 

Raw HiFi reads were downloaded from the NCBI Sequence Read Archive for *Cebidichthys violaceus* SRR19653852 (Wright et al., 2023); *Cryptacanthodes maculatus* SRR29012721; *Leptoclinus maculatus* SRR25603844, SRR25603845; *Melanostigma gelatinosum* ERR10879927 (Bista et al., 2024); and *Pholis gunnellus* ERR6436364, ERR6412365 (Programme et al., 2022). For *Lycodopsis pacificus*, FASTQ files containing raw HiFi reads were retrieved from GenomeARK. Hi-C data for *C. maculatus* (SRX26356541) was also downloaded from the from the NCBI Sequence Read Archive. 

Raw HiFi BAM and FASTQ for *Lycodicthys dearborni* and *Zoarces americanus* are available on NCBI SRA under the accession PRJNA1236397 [to be unembargoed upon acceptance]. HiFi FASTQ reads can be accessed on SRA for *Anarhichas lupus* (PRJNA980960) and *Anarhichas minor* (PRJNA982125) [to be unembargoed upon acceptance]

## Table of contents

### 1. genome_assembly_annotation
**sra_accession.sh**: Bash script used to download raw reads and Hi-C data from NCBI Sequence Read Archive using sratoolkit v3.0.0. 

**hifiasm.sh**: Bash script for assembly of raw Hifi reads using hifiasm v0.19.9. 

**hifiasm_hi_C.sh**: Bash script for assembly of *Cryptacanthodes maculatus* genome with Hi-C data using hifiasm v0.19.9. 

**trimgalore_hi_c.sh**: Bash script used to trim raw *Cryptacanthodes maculatus* Hi-C data before assembly using trimgalore v0.6.10

**AFP_annotation.sh**: Bash script for conversion of hifiasm output .gfa to .fasta and annotation of AFP III genes. haplotype 1 and haplotype 2 phased .gfa files and the raw unitig .gfa are converted to .fasta format using gfatools v0.5. Then, annotations of AFP III genes are performed using Exonerate v2.4.0 with a query sequence of the translated *Macrozoarces americanus* AFP (Mamericanus_AFP.txt)

**flanking_gene_annotation.sh**: Bash script for the annotation of flanking genes to the ancestral AFP array in Zoarcoidei. Annotation of the upstream *sncgb* gene in the haplotype-phased and unitig .fasta is performed using Exonerate v2.4.0 with a query *Gasterosteus aculeatus sncgb* protein sequence (Gacul_sncgb.fa). The same annotation is then performed with a query protein sequence of the *G. aculetaus ldg3b* gene (Gacul_ldb3b.fa), which is downstream of the ancestral array. 

### 2. analysis
**switch_error_screen.sh**:



### 3. figure_code
**haplo_afp.Rmd**: R markdown file containing all code used for figures 5 (species tree of Zoarcoidei and dotplot of uncertatinty in AFP copy number created by misassembly) and S1 (correlation between median copy number for each AFP array and uncertainty in copy number) using copy number counts calculated by *gfa_parser* (zoarcoidei_haplo-afp - fig_5.csv). 

