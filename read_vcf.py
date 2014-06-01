import numpy as np
import vcf
import matplotlib.pyplot as plt

vcf_reader= vcf.Reader(open("act.raw.vcf","r"))
vcf_cont=[]
for rec in vcf_reader:
    vcf_cont.append(rec)
vcf_10=[]
for rec in vcf_cont:
    i=0
    for sample in rec.samples:
        if sample.data[0]==None:
	    next
        elif sum(sample.data[1])>=10:
	    i+=1
    if i==len(vcf_cont[0].samples):
        vcf_10.append(rec)
writer=vcf.Writer(open("act.10.vcf","w"), vcf.Reader(open("act.raw.vcf","r")))
for i in range(0, len(vcf_10)):
    writer.write_record(vcf_10[i])
writer.close()
