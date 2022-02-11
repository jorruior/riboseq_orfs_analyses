##Requirements: PRANK
#Argument 1: Output folder of script 1, including all maf alignments
#Argument 2: Tree file, in nwk format
#Argument 3: Output tag

for f in $1/orfs/*.maf; do echo $(basename $f); prank -d=$f -showanc -showevents -prunetree -F -once -t=$2 -o=tmp/$(basename $f); done

rm $3.ancestors
echo -e 'orf_id\tsp\tev_age\tsyn_age\tgained\tgained_convergent\tlost\tdenovo\tseq' > $3.ancestors
for f in tmp/*best.anc.dnd; do python3 parsing_ancestors.py $f $3; done
