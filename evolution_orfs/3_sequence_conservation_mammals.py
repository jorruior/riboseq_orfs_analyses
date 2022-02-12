import sys
import string
import subprocess
import os
import glob
import statistics
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from optparse import OptionParser

__author__ = "Jorge Ruiz-Orera"
__contributor__="..."
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Jorge Ruiz-Orera"
__email__ = "jorruior@gmail.com"

#Script to aligned ORF against orthologous regions. Requirements: BLAST

def check_arg (arg_str,s):
	'''Check if arg was written'''
	if not arg_str:   # if filename is not given
		print("Error: " + str(s) + " argument not given\n")
		exit()

#Arguments
usage = "\n%prog [options]"
parser = OptionParser(usage,version="%prog " + __version__)
parser.add_option("-i","--input",action="store",dest="maf_folder",help="(Required) Folder with all MAF alignments (output of 1_extract_multiple_alignments.py)")
parser.add_option("-p","--prot",action="store",dest="prot",help="(Required) Fasta including all translated sequences encoded by the ORFs")

(opt,args)=parser.parse_args()

check_arg(opt.maf_folder,"--maf_folder")
check_arg(opt.maf_folder,"--prot")

#Main
primatomorpha = ("hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1","colAng1","HLpilTep1","calJac3","aotNan1","saiBol1","cebCap1","tarSyr2","otoGar3","micMur3","proCoq1", "galVar1","tupChi1","eulMac1","eulFla1")

total_sp = []
folder = opt.maf_folder + "/orfs/*.maf"
for file in glob.glob(folder):
		seqs = SeqIO.index(file, "fasta")
		for species in seqs:
			if species == "hg38":
				continue
			if not species in total_sp:
				total_sp.append(species)

os.system("mkdir " + opt.maf_folder + "/fastas")
outs = {}
for sp in total_sp:
	outs[sp] = open(opt.maf_folder + "/fastas/" + sp + "_aligned_proteins.fa","w+")

folder = opt.maf_folder + "/orfs/*.fa"
for file in glob.glob(folder):
	for line in open(file):
			if ">" in line:
				s = line.split()[0].replace(">","")
				species = line.split()[2]
				if species != "hg38":
					outs[species].write(">" + s + "\n")
			else:
				if species == "hg38":
					exc = []
					for n,c in enumerate(line.rstrip("\n")[:-1]):
						if c == "*":
							exc.append(n)
					continue

				outs[species].write(line.rstrip("\n") + "\n")
				trunc = 0
				l = 0
				for n,c in enumerate(line.rstrip("\n")[:-1]):
					if not n in exc:
						if c == "*":
							trunc += 1
						l += 1

cs = {}
for sp in total_sp:
	if sp in primatomorpha: #Primatomorpha is not evaluated
		continue
	outs[sp].close()
	print("BLAST to " + sp)
	os.system("blastp -query " + opt.prot + " -subject " + opt.maf_folder + "/fastas/" + sp + "_aligned_proteins.fa -evalue 1e-2 -max_target_seqs 1 -outfmt \"6 qseqid sseqid nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" >  " + opt.maf_folder + "/fastas/blast_all_orfs_to_" + sp + ".out")
	for line in open(opt.maf_folder + "/fastas/blast_all_orfs_to_" + sp + ".out"):
		name = line.split("\t")[0]
		if not name in cs:
			cs[name] = [[],[]]
		evalue = float(line.split("\t")[10])
		cs[name][0].append(sp)
		cs[name][1].append(evalue)


out = open(opt.maf_folder + "/fastas/conservation_scores.tsv","w+")
out.write("orf_id\tCS\tall_evalues\tall_species\n")
for name in cs: #Only cases with at least one match are reported
	out.write(name + "\t" + str(-math.log(statistics.median(cs[name][1]),10)) + "\t" + ",".join(map(str,cs[name][1])) + "\t" + ",".join(map(str,cs[name][0])) + "\n")
out.close()

exit()
