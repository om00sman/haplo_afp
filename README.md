# haplo_afp

## Authors

Owen Moosman, University of California, Santa Cruz, Department of Molecular, Cell, & Developmental Biology

Samuel Bogan, University of California, Santa Cruz, Department of Ecology and Evolutionary Biology

Joanna Kelley, University of California, Santa Cruz, Department of Ecology and Evolutionary Biology

## Description

This repository contains code used by Owen Moosman, Samuel Bogan, and Joanna Kelley for analyzing copy number variation of antifreeze protein (AFP) genes between haplotypes of polar fish in the suborder Zoarcoidei. Bioinformatic tools developed and used for this analysis can be found elsewhere. *gfa_parser* is hosted and described at the Github repository https://github.com/snbogan/gfa_parser. *switch_error_screen* is hosted and described at the Github repository https://github.com/snbogan/switch_error_screen. All R scripts were run using R version 4.3.2. 

This research was funded by the United States National Science Foundation, Office of Polar Programs (Award number 2312253) 

## Data Availabilty 

Raw HiFi reads were downloaded from the NCBI Sequence Read Archive for *Leptoclinus maculatus* SRR25603844, SRR25603845; and *Pholis gunnellus* ERR6436364, ERR6412365 (Programme et al., 2022). Raw HiFi reads for the wolffishes Anarhichas lupus and Anarhichas minor were obtained from Bogan et al. (Bogan et al. 2024).

Raw HiFi BAM and FASTQ for *Lycodicthys dearborni* and *Zoarces americanus* are available on NCBI SRA under the accession PRJNA1236397 [to be unembargoed upon acceptance]. HiFi FASTQ reads can be accessed on SRA for *Anarhichas lupus* (PRJNA980960) and *Anarhichas minor* (PRJNA982125) [to be unembargoed upon acceptance]

## Table of contents

### 1. genome_assembly_annotation
**sra_accession.sh**: Bash script used to download raw reads and Hi-C data from NCBI Sequence Read Archive using sratoolkit v3.0.0. 

**hifiasm.sh**: Bash script for assembly of raw Hifi reads using hifiasm v0.19.9. 

**afp_annotation.sh**: Bash script for conversion of hifiasm output .gfa to .fasta and annotation of AFP III genes. Haplotype 1 and haplotype 2 phased .gfa files and the raw unitig .gfa are converted to .fasta format using gfatools v0.5. Then, annotations of AFP III genes are performed using Exonerate v2.4.0 with a query sequence of the translated *Zoarces americanus* AFP (Mamericanus_AFP.txt)

**afp_consensus_compare.txt**: Contains code to generate a consensus AFP sequence for each species based on exonerate hits generated with the *Zoarces americanus* americanus query and perform exonerate annotations with this consensus sequence. 

**flanking_gene_annotation.sh**: Bash script for the annotation of flanking genes to the ancestral AFP array in Zoarcoidei. Annotation of the upstream *sncgb* gene in the haplotype-phased and unitig .fasta is performed using Exonerate v2.4.0 with a query *Gasterosteus aculeatus sncgb* protein sequence (Gacul_sncgb.fa). The same annotation is then performed with a query protein sequence of the *G. aculetaus ldg3b* gene (Gacul_ldb3b.fa), which is downstream of the ancestral array. 

### 2. analysis
**switch_error_screen.sh**: A bash script that extracts read phasing information from the output .gfa's of haplotypes 1 and 2, then performs a realignment of all reads to a combined assembly of haplotypes 1 and 2 using pbmm2 v1.14.99 (https://github.com/PacificBiosciences/pbmm2), a wrapper of minimap2 optimized for alignment of HiFi reads. The resulting .bam is filtered into two .bam's using samtools v1.21, one containing reads phased to haplotype 1 and unphased reads, and the other containing reads phased to haplotype 2 and unphased reads. Finally, an Exonerate query (as described above) is used to annotate the combined haplotype 1 and 2 assembly for use in alignment visualization software. 

**AFP_filter.sh**: A bash script that extracts contigs containing AFPs (as determined by exonerate annotation) from the larger .bam file using samtools v1.21 (in order to reduce file size). 

**dearborni_switch_err.sh**: A bash script similar to **switch_error_screen.sh** but that calculates the percent of reads mapped to the combined assembly and also performs a less stringent mapping. 

**pgunn_align.sh**: A bash script that performs an alignment of the raw haplotype 1 *Pholis gunnellus* reads to a reference with a simulated switch error. 

### 3. figure_code
**haplo_afp.Rmd**: R markdown file containing all code used for figures 6 (species tree of Zoarcoidei and dotplot of uncertatinty in AFP copy number created by misassembly) and S1 (correlation between median copy number for each AFP array and uncertainty in copy number) using copy number counts calculated by *gfa_parser* (zoarcoidei_haplo-afp - fig_5.csv). 

