# Analyzing the evolution of ORF sequences and structures across primate and mammalian species
This repository contains a series of scripts developed to obtain the main results of the manuscripts "Origins, mechanisms, and implications of de novo gene and ORF birth in humans" and "Evolutionary origins and interactomes of human young microproteins and small peptides translated from short open reading frames".

**-evolution_orfs:** Scripts to build multiple alignments of ORF sequences, reconstruct ancestral sequences, and evaluate the sequence conservation of ORFs across the included species. 
```
bash pipeline.sh [BED ORFs 1-based] [folder with MAF alignments: chr1.maf, chr2.maf...] [output folder] [nwk file] [FASTA ORF proteins]
```
