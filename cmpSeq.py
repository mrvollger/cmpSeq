import os
import glob
from Bio import SeqIO
import re
import pandas as pd
import itertools
pd.set_option('display.width', 1000)


SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
#shell.prefix("source %s/env_PSV.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)
#shell.prefix("source %s/env_PSV.cfg; " % SNAKEMAKE_DIR)
cwd=os.getcwd()
LAdir = "/".join(os.getcwd().split("/")[:-1])
base2="/net/eichler/vol2/home/mvollger/projects/crossMatch"


# setup pairwise widcards 
pair1=[]
pair2=[] 
pairs = []
with open("cmp.fasta", "r") as handle:
	records = list(SeqIO.parse(handle, "fasta"))
	for pair in itertools.combinations(records, 2):
		name = pair[0].id + "__" + pair[1].id
		pair1.append(pair[0].id)
		pair2.append(pair[1].id)
		pairs.append(name)
print(pairs)


rule all:	
	input:
		plot="alignment.png",
	shell:
		"""
		echo done	
		"""

#
# split fastas 
#
rule split:
	input:
		fasta = "cmp.fasta",
	output:
		fasta = expand("pairfasta/{pair}.fasta", pair=pairs)
	run:	
		shell("mkdir -p pairfasta")
		with open(input["fasta"], "r") as handle:
			records = list(SeqIO.parse(handle, "fasta"))
			for pair in itertools.combinations(records, 2):
				name = pair[0].id + "__" + pair[1].id
				SeqIO.write(pair, "pairfasta/{}.fasta".format(name), "fasta")

#
# miropeats 
#
rule miropeats:
	input:
		fasta = "pairfasta/{pair}.fasta",
	output:
		miro="miropeats/{pair}.pdf",
	params:
		threshold="500",
	shell: 
		"""
		mkdir -p miropeats
		cd miropeats
		miropeats -s {params.threshold} -onlyinter ../{input.fasta} > contig_compare
		#miropeats -s {params.threshold} ../{input.fasta} > contig_compare
		if [ -f threshold{params.threshold} ]; then
			mv threshold{params.threshold} {wildcards.pair}.ps
			ps2pdf {wildcards.pair}.ps
		else
			touch ../{output.miro}
		fi 

		"""



#
# alignments 
#
rule alignments:
	input:
		miro="miropeats/{pair}.pdf",
		fasta = "pairfasta/{pair}.fasta",
	output:
		aln="align/{pair}.aln",
	threads: 8
	shell:
		"""
		mkdir -p align 
		export PATH=$PATH:/net/eichler/vol2/local/inhousebin
		cross_match -alignments -masklevel 0 -minmatch 100 -tags \
				-bandwidth 500 -gap_ext -1 -gap_init -10 -penalty -1 \
				{input.fasta} > {output.aln}
		"""


#
# new perID  
#
rule perID:
	input:
		alnsam=expand("align/{pair}.aln", pair=pairs),
	output:
		allaln="align/all.aln",
		alnsam="align/aln.sam",
		perid="perID/slidePerID.tsv",
		plot="alignment.png",
	shell:
		"""
		mkdir -p perID
		cat {input.alnsam} > {output.allaln}
		cmToSam.py {output.allaln} {output.alnsam} 
		samPerID.py --header --window 1000 --step 200 {output.alnsam} {output.perid}
		slideAln.R {output.perid} {output.plot}
		"""	









