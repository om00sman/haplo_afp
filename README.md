Authors
---

Owen Moosman, University of California, Santa Cruz, Department of Molecular, Cell, & Developmental Biology

Samuel Bogan, University of California, Santa Cruz, Department of Ecology and Evolutionary Biology

Joanna Kelley, University of California, Santa Cruz, Department of Ecology and Evolutionary Biology

Description
---

This repository contains code used by Owen Moosman, Samuel Bogan, and Joanna Kelley for analyzing copy number variation of antifreeze protein genes between haplotypes of polar fish in the suborder Zoarcoidei. Bioinformatic tools developed and used for this analysis can be found elsewhere. gfa_parser is hosted and described at the Github repository https://github.com/snbogan/gfa_parser. switch_error_screen is hosted and described at the Github repository https://github.com/snbogan/switch_error_screen.

This research was funded by the United States National Science Foundation, Office of Polar programs award number 2312253 

Data Availabilty 
---

Raw HiFi reads were downloaded from the NCBI Sequence Read Archive for Cebidichthys violaceus SRR19653852 (Wright et al., 2023); Cryptacanthodes maculatus SRR29012721; Leptoclinus maculatus SRR25603844, SRR25603845; Melanostigma gelatinosum ERR10879927 (Bista et al., 2024); and Pholis gunnellus ERR6436364, ERR6412365 (Programme et al., 2022). For Lycodopsis pacificus, FASTQ files containing raw HiFi reads were retrieved from GenomeARK. Hi-C data for C. maculatus (SRX26356541) was also downloaded from the from the NCBI Sequence Read Archive. 

Raw HiFi BAM and FASTQ for Lycodicthys dearborni and Zoarces americanus are available on NCBI SRA under the accession PRJNA1236397 [to be unembargoed upon acceptance]. HiFi FASTQ reads can be accessed on SRA for Anarhichas lupus (PRJNA980960) and Anarhichas minor (PRJNA982125) [to be unembargoed upon acceptance]

Table of contents
---

1. genome_assembly_annotation
     i. sra_accession.sh: Bash script used to download raw reads and Hi-C data from NCBI Sequence Read Archive. 
3. analysis
4. figure_code
