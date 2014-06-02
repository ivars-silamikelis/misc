!#/usr/bin/python
# -*- coding: utf-8 -*-
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from Bio import SeqIO
coords={"starts":[],"ends":[]}
for path in glob.glob("CTRI-2.Russia.CP002992/*.crd"):
	name_and_ref=os.path.splitext(os.path.basename(path))[0]
	ref=os.path.dirname(path).split("/")[-1]    #ja izmanto absolute path
	name=name_and_ref.replace("."+ref,"")
	ends=get_overlaps(path)
	coords["starts"]=coords["starts"]+ends["starts"]
	coords["ends"]=coords["ends"]+ends["ends"]

for seq in SeqIO.parse("/home/ivars/Data/Mtb/references/no_IS/CTRI-2.Russia.CP002992_clean.fasta","fasta"):
	#print seq.id
tsds=[]
for i in range(0, len(coords["starts"])):
	tsds.append( seq.seq.tostring()[coords["starts"][i]:coords["ends"][i]])
max_len=np.max([len(x) for x in tsds])

string_stats=dict()
for i in range(0, max_len):
	string_stats[i]={"A":0,"T":0,"G":0,"C":0}

for tsd in tsds:
	for i in range(0,len(tsd)):
		string_stats[i][tsd[i]]+=1
		print string_stats
def get_overlaps(crd_file):
	coords=dict()
	coords["starts"]=[]
	coords["ends"]=[]
	with open(crd_file) as f:
		for i, l in enumerate(f):
			other_start=int(l.split("\t")[0])
			start=int(l.split("\t")[1])
			end=int(l.split("\t")[2])
			other_end=int(l.split("\t")[3])
			if start > end and other_start < end and other_end > start and other_start < other_end:
				coords["starts"].append(end)
				coords["ends"].append(start)
			if start < end and other_start > end and other_end < start and other_start > other_end:
				coords["starts"].append(start)
				coords["ends"].append(end)
	return coords

