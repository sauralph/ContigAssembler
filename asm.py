#!/bin/python
# Contig Assembler
from human_format import human_format
from Bio import SeqIO
import os
import random
import pandas as pd



#DEBUG
ignore = 0

contigs    = SeqIO.parse("v80.fasta", "fasta")
references = SeqIO.parse("ref.fasta", "fasta") 
reference  = next(references)
#What is the reference?
print("Using as reference {0}".format(reference.id))
ref_size   = len(reference.seq)
print("Reference Size: {0}".format(human_format(ref_size)))
#Foreach Contig
for contig in contigs:
	ignore = ignore +1
	if(ignore <7 ):
		continue

	#Show Contig Name
	print("Contig {0}".format(contig.name))
	#Show Contig Size
	print("Size: {0}pb".format(human_format(len(contig.seq))))
	#create a Fasta file for current contig
	SeqIO.write(contig,"current.fasta","fasta")
	run_scores = []
	run_contig = []
	run_end    = []
	run_start  = []
	#Sample Ref Genome And Align
	for i in range(1,100):
		start_base = random.randint(1,(ref_size-500))
		end_base   = start_base + 500
		outfile    = "{0}-{1}-{2}".format(contig.name,start_base,end_base)
		#execute Waterman Algorithm
		os.system("water -asequence current.fasta -bsequence ref.fasta -gapopen 10 -gapextend 0.5 -sbegin2 {0} -send2 {1} -outfile {2}".format(
			start_base,
			end_base,
			outfile
			))
		score = os.popen("grep Score {0}|cut -f3 -d' '".format(outfile)).read()
		score = score.splitlines()
		score = score[0]
		#print(score)
		run_scores.append(score)
		run_contig.append(contig.name)
		run_end.append(end_base)
		run_start.append(start_base)
	df = pd.DataFrame({
		"contig":run_contig,
		"scores":run_scores,
		"end":run_end,
		"start":run_start
	})
	#Print Start and End of HIT zone and save it to CSV
	df.to_csv("contig{0}.csv".format(ignore))