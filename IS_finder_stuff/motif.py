#!/usr/bin/python
# -*- coding: utf-8 -*-
import glob
#import matplotlib.pyplot as plt
#import numpy as np
#import pandas as pd
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
	print seq.id
tsds=[]
for i in range(0, len(coords["starts"])):
	tsds.append( seq.seq.tostring()[coords["starts"][i]:coords["ends"][i]])
max_len=np.max([len(x) for x in tsds])

string_stats=dict()
for i in range(0, max_len):
	string_stats[i]={"A":0,"T":0,"G":0,"C":0}

counts=collections.Counter([len(x) for x in tsds])

for tsd in tsds:
	for i in range(0,len(tsd)):
		string_stats[i][tsd[i]]+=1
		#print string_stats

ref_stats=get_string_stats(sequence)
summa=ref_stats["A"]+ref_stats["T"]+ref_stats["C"]+ref_stats["G"]
for nuc in ref_stats.keys():
    print nuc+":"+str(ref_stats[nuc]/float(summa))

for i in range(0,max_len):
	summa=string_stats[i]["A"]+string_stats[i]["T"]+string_stats[i]["C"]+string_stats[i]["G"]
	for nuc in string_stats[i].keys():
		print nuc+str(i)+":"+str(string_stats[i][nuc]/float(summa))



def get_string_stats(sequence):
	string_stats={"A":0,"T":0,"G":0,"C":0}
	for i in sequence:
		string_stats[i]+=1
	return string_stats
	
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

