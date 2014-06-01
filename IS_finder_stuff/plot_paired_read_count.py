#!/usr/bin/python
# -*- coding: utf-8 -*-
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import matplotlib.cm as cm
data=dict()
paths=sorted(glob.glob("*/S*.crd"))
#inicializē mainīgo
for path in paths:
	name_and_ref=os.path.splitext(os.path.basename(path))[0]
	ref=os.path.dirname(path).split("/")[-1]    #ja izmanto absolute path
	name=name_and_ref.replace("."+ref,"")
	data[name]=dict()
	#data[name]["count"]=[]
	data[name]=[]
	#data[name]["ref"]=[]
data["ref"]=[]
#Datu iegūšana
directories=glob.glob("*/")

for directory in directories:
	ref=os.path.dirname(directory).split("/")[-1]
	data["ref"].append(ref)
	paths=glob.glob(ref+"/S*.crd")
	for path in paths:
		#path=names[0]
		name_and_ref=os.path.splitext(os.path.basename(path))[0]
		ref=os.path.dirname(path).split("/")[-1]	#ja izmanto absolute path
		name=name_and_ref.replace("."+ref,"")
		#data[name]["count"].append(line_count(path))
		#data[name]["ref"].append(ref)
		data[name].append(line_count(path))
df=pd.DataFrame(data)
krasas=["","#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5"]#8dd3c7
#ffffb3
#bebada
#fb8072
#80b1d3
#fdb462
#b3de69
#fccde5
bar_width = 0.35
i=1
index=np.arange(0,len(df["ref"]))*bar_width*(len(df.keys())+1)
for sample in df.keys():
	if sample=="ref":
		continue
	plt.bar(index+bar_width*i,df[sample],bar_width, label=sample, color=krasas[i])#cm.hot(i/10.,1))
	i+=1

plt.xticks(index + bar_width, df["ref"])
plt.legend()

plt.tight_layout()
plt.show()


def line_count(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
	    try:
			i
		except NameError:
			return 0
		else:
			return i + 1
