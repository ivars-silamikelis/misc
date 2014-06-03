import vcf
import numpy as np
import pandas as pd
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
def read_vcf(name):
    vcf_raw=[]
    vcf_reader= vcf.Reader(open(name,"r"))
    for rec in vcf_reader:
        vcf_raw.append(rec)
    return vcf_raw

sequences=dict()
seqs_string=dict()
vcf_content=read_vcf("output.vcf")
#inicialize dictionary
for sample in vcf_content[0].samples:
	sam=re.sub(".variant\d*","",sample.sample)
	sequences[sam]=[]

for rec in vcf_content:
	for sample_rec in rec.samples:
		sam=re.sub(".variant\d*","",sample_rec.sample)
		if sample_rec["GT"]==None or int(sample_rec["GT"])==0:
			sequences[sam].append(rec.REF)
		else:
			gt=int(sample_rec["GT"])-1
			sequences[sam].append(str(rec.ALT[gt]))
for sample in sequences.keys():
	seqs_string[sample]="".join(sequences[sample])
records=[]
for sample in seqs_string:
	records.append(SeqRecord(Seq(seqs_string[sample]), id=sample, description=""))
SeqIO.write(records, "output.fasta", "fasta")

