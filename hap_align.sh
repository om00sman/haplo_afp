
#load required modules

module load minimap2

minimap2 -ax asm5 ref.fa asm.fa > aln.sam
