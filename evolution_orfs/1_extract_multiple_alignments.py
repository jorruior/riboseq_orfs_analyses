#!/usr/bin/env python
import sys
import string
import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.AlignIO import MafIO
from optparse import OptionParser

__author__ = "Jorge Ruiz-Orera"
__contributor__="..."
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Jorge Ruiz-Orera"
__email__ = "jorruior@gmail.com"

#Script to get all ORF alignments based on MAF file. Requirements: Biopython, PHAST (maf_parse)

class trans_object:
	def __init__(self, chrm, gene, strand, start, end):
		self.chrm = chrm
		self.gene = gene
		self.strand = strand
		self.start = start
		self.end = end


def check_arg (arg_str,s):
	'''Check if arg was written'''
	if not arg_str:   # if filename is not given
		print("Error: " + str(s) + " argument not given\n")
		exit()


def check_file (file_str):
	'''Check if input really exists'''
	try:
		open("%s" %file_str)
	except:
		print("Error: " + file_str + " input not found\n")
		exit()


def parse_bed(bed):
	#Read a bed and create a dict with sorted transcript coordinates, chrm, strand, and gene
	trans = {}
	for line in open(bed):
		t_name = line.split('\t')[3]
		chrm = line.split("\t")[0]
		if not "chr" in chrm:
			chrm = "chr" + chrm
		
		trans.setdefault(t_name,trans_object(chrm,t_name,line.split("\t")[5].rstrip('\n'),[],[]))
		trans[t_name].start.append(int(line.split("\t")[1]))
		trans[t_name].end.append(int(line.split("\t")[2]))

	[trans[x].start.sort() for x in trans]
	[trans[x].end.sort() for x in trans]	

	return trans


def extract_sequences(ma,pin,pif,regin,focal):
	all_seqs = {}
	for s in ma:
		if not s.name.split(".")[0] in all_seqs:
			all_seqs[s.name.split(".")[0]] = ["",0,0,s.name.split(".")[1],s.annotations["start"],s.annotations["strand"]]
		if s.name.split(".")[0] == focal:
			fseq = ""
			counter = regin
			p1 = -1
			p2 = -1
			p11 = 0
			p22 = 0
			for n,ss in enumerate(str(s.seq)):
				if (counter >= pin - 1) and (counter <= pif - 1):
					if p1 == -1:
						p1 = n
					fseq = fseq + ss
				elif (counter > pif - 1) and (p2 == -1):
					p2 = n-1
				if ss != "-":
					counter += 1
			if p2 == -1:
				p2 = n
		else:
			fseq = ""
			p11 = 0
			p22 = 0
			for n,ss in enumerate(str(s.seq)):
				if (n < p1) and (ss != "-"):
					 p11 += 1
				if (n < p2) and (ss != "-"):
					 p22 += 1
				if (n >= p1) and (n <= p2):
					fseq = fseq + ss				

		all_seqs[s.name.split(".")[0]] = [fseq,p11,p22,s.name.split(".")[1],s.annotations["start"],s.annotations["strand"]]

	return all_seqs


def get_seq_maf(orf,maf,focal,tag):
	#It uses maf_parse from the PHAST utilities
	orf_seq = {}
	coords = []
	l = 0
	for n,sto in enumerate(orf.start):	
		eno = orf.end[n]
		m = "none"
		op = 0
		l = l + orf.end[n] - orf.start[n] + 1
		os.system("maf_parse " + maf + " -s " + str(orf.start[n]) + " -e " + str(orf.end[n]) + " -o MAF > tmp/" + tag + "-" + orf.gene + ".subset.maf")
		for ma in AlignIO.parse("tmp/" + tag + "-" + orf.gene + ".subset.maf", "maf"):
			st = int(ma[0].annotations["start"])
			strand = ma[0].annotations["strand"]
			usize = int(ma[0].annotations["size"]) 
			tsize = int(ma[0].annotations["srcSize"])
			if str(st) + "/" + str(usize) in coords:
				continue
			if (sto >= st) and (sto <= st+usize):
				op = 1
				all_seqs = extract_sequences(ma,sto,eno,st,"hg38")
				coords.append(str(st) + "/" + str(usize))
				
				for sp in all_seqs:
					if not sp in orf_seq:
						orf_seq[sp] = ["",[],[],[],[]]
					orf_seq[sp][0] = orf_seq[sp][0] + all_seqs[sp][0]
					orf_seq[sp][1].append(all_seqs[sp][4] + all_seqs[sp][1])
					orf_seq[sp][2].append(all_seqs[sp][4] + all_seqs[sp][2])
					orf_seq[sp][3].append(all_seqs[sp][3])
					orf_seq[sp][4].append(all_seqs[sp][5])

				if (eno >= st) and (eno <= st+usize):
					break
			elif (op == 1) and (eno >= st) and (eno <= st+usize):
				all_seqs = extract_sequences(ma,sto,eno,st,"hg38")
				coords.append(str(st) + "/" + str(usize))

				for sp in all_seqs:
					if not sp in orf_seq:
						orf_seq[sp] = ["",[],[],[],[]]
					orf_seq[sp][0] = orf_seq[sp][0] + all_seqs[sp][0]
					orf_seq[sp][1].append(all_seqs[sp][4] + all_seqs[sp][1])
					orf_seq[sp][2].append(all_seqs[sp][4] + all_seqs[sp][2])
					orf_seq[sp][3].append(all_seqs[sp][3])
					orf_seq[sp][4].append(all_seqs[sp][5])

				break
			elif (op == 1) and (eno >= st) and (eno > st+usize):
				coords.append(str(st) + "/" + str(usize))
				all_seqs = extract_sequences(ma,sto,eno,st,"hg38")

				for sp in all_seqs:
					if not sp in orf_seq:
						orf_seq[sp] = ["",[],[],[],[]]
					orf_seq[sp][0] = orf_seq[sp][0] + all_seqs[sp][0]
					orf_seq[sp][1].append(all_seqs[sp][4] + all_seqs[sp][1])
					orf_seq[sp][2].append(all_seqs[sp][4] + all_seqs[sp][2])
					orf_seq[sp][3].append(all_seqs[sp][3])
					orf_seq[sp][4].append(all_seqs[sp][5])

			elif (eno < st):
				break
		os.system("rm tmp/" + tag + "-" + orf.gene + ".subset.maf")

	return orf_seq,l,coords

#Arguments
usage = "\n%prog, we recommend to parallelize this script in a computing cluster since maf alignments are slow to parse [options]"
parser = OptionParser(usage,version="%prog " + __version__)
parser.add_option("-b","--bed",action="store",dest="bed",help="(Required) 1-based BED file with ORF strand and coordinates")
parser.add_option("-m","--maf_folder",action="store",dest="maf_folder",help="(Required) Folder storing all maf files separated by chromosome (chr1.maf, chr2.maf...)")
parser.add_option("-o","--output",action="store",dest="tag",help="(Required) Name of output folder")
parser.add_option("-f","--force",action="store",dest="force",default="no",help="If 'yes', force to re-run alignments that are already calculated in the same folder. (default = 'no')")

(opt,args)=parser.parse_args()

check_arg(opt.bed,"--bed")
check_arg(opt.maf_folder,"--maf_folder")
check_arg(opt.tag,"--output")
check_file(opt.bed)

try:
	os.system("mkdir " + opt.tag)
	os.system("mkdir " + opt.tag + "/orfs")
	os.system("mkdir " + opt.tag + "/beds")
except:
	print("No permissions to write folder, exiting.")
	exit(0)

#Main
orfs = parse_bed(opt.bed)
for count,orf in enumerate(orfs):
	if count % 5 == 0:
		print(str(count) + " out of " + str(len(orfs)))
	if opt.force == "yes":
		if os.path.isfile(opt.tag + "/orfs/" + orf.replace(":","__") + ".maf"):
			print (orf + " file exist")
			continue
	(orf_seq,l,coords) = get_seq_maf(orfs[orf],opt.maf_folder + "/" + orfs[orf].chrm + ".maf","hg38", opt.tag)
	if not "hg38" in orf_seq:
		print(orf + "\t0\n", end = "")
		continue
	m = str(float(len(orf_seq["hg38"][0].replace("-",""))/float(l)))
	out1 = open(opt.tag + "/orfs/" + orf.replace(":","__") + ".fa","w+")
	out2 = open(opt.tag + "/orfs/" + orf.replace(":","__") + ".maf","w+")
	out3 = open(opt.tag + "/beds/" + orf.replace(":","__") + ".bed","w+")
	for sp in orf_seq:
		sq = orf_seq[sp][0]
		if orfs[orf].strand == "-":
			sq = str(Seq(orf_seq[sp][0]).reverse_complement())
		if len(sq) == len(orf_seq["hg38"][0]):
			out1.write(">" + orf + " " + orfs[orf].gene  + " " + sp + " " + m + "\n" + str(Seq(sq.replace("-","")).translate(cds=False)) + "\n")
			out2.write(">" + sp + "\n" + sq + "\n")
			t2 = ""
			if(len(set(orf_seq[sp][3]))!=1):
				t2 = "_splitted"
			if sp != "hg38":
				for n,elemento in enumerate(orf_seq[sp][1]):
					out3.write(orf_seq[sp][3][n] + "\t" + str(orf_seq[sp][1][n]) + "\t" + str(orf_seq[sp][2][n]) + "\t" + orf + "\t" + sp + t2 + "\t" + str(orf_seq[sp][4][n]).replace("-1","-").replace("1","+") + "\n")

	out1.close()
	out2.close()
	out3.close()

exit(0)
