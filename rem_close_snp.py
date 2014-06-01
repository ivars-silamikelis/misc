#!/usr/bin/python
#import numpy as np
#import pandas as pd
vcfs=["mtb_h37rv.raw"]
thresh=50
import re
for sn in vcfs:
    vcf=[]
    hdr=[]
    with open(sn+".vcf") as infile:
        for line in infile:
            if not re.search('^#',line):
    	        vcf.append((re.sub('\n','',line)).split("\t"))
            else:
	        hdr.append(line)

    pos=[int(x[1]) for x in vcf]
    good_idx=[]
    if (pos[1]-pos[0])>thresh:
        good_idx.append(0)

    for i in range(1, len(pos)-1):
        if (pos[i]-pos[i-1]>thresh) and (pos[i+1]-pos[i]>thresh):
            good_idx.append(i)
    if (pos[-1]-pos[-2])>thresh:
        good_idx.append(len(pos)-1)
   
    good_vcf=["\t".join(vcf[x])+"\n" for x in good_idx]
    final=hdr+good_vcf
    f=open(sn+".near"+str(thresh)+".vcf","w")
    for x in final:
        f.write(x)

    f.close()


