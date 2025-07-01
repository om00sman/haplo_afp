species=l_maculatus
in=/hb/groups/kelley_training/owen/zoarcoidei/data/assemblies/$species
out=/hb/home/omoosman/owen/zoarcoidei/analysis/complexity

gfatools view -l $(cat a_minor_trans1.txt) $in/${species}_hifi.asm.bp.r_utg.gfa > a_minor_trans1.gfa
sed 's/, /,/g' a_minor_array.txt > a_minor_trans1_.txt
