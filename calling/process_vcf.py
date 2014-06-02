import vcf
import numpy as np
import pandas as pd
import glob

reg1=pd.read_table("rep_ppe_reg.txt")
reg2=pd.read_table("all.txt")
cover_thresh=10
saraksts=glob.glob("*.vcf")
names_pop=[x.replace(".samtools.vcf","").replace(".gatk_ug.vcf","").replace(".varscan.vcf","") for x in saraksts]
names_uniq=set(names_pop)
for name in names_uniq:
	print name
	sam=read_vcf(name+".samtools.vcf")
	gatk=read_vcf(name+".gatk_ug.vcf")
	varscan=read_vcf(name+".varscan.vcf")
	vs_homs=[x for x in varscan if x.INFO["HOM"]==1]
	#filtrs pec coverage un homozigotates
	gatk_homs=info_cover_filter(gatk,cover_thresh)
	sam_homs=info_cover_filter([x for x in sam if x.samples[0]["GT"]=="1/1"],cover_thresh)
	varscan_homs=cover_filter([x for x in varscan if x.INFO["HOM"]==1], cover_thresh)
	gatk_final=keep_regs(filter_regs(gatk_homs, reg1),reg2)
	sam_final=keep_regs(filter_regs(gatk_homs, reg1), reg2)
	varscan_final=keep_regs(filter_regs(vs_homs, reg1),reg2)

def filter_regs(snps, regs):
	good_snps=snps
	is_in_region=0	
	for snp in good_snps:
		for i in range(0, len(regs["start"])):
			if snp.POS >= regs["start"][i] and snp.POS<=regs["end"][i]:
				is_in_region=1
		if is_in_region==1:
			good_snps.remove(snp)
			is_in_region=0
	return good_snps
def keep_regs(snps, regs):
	good_snps=[]
	is_in_region=0	
	for snp in snps:
		for i in range(0, len(regs["start"])):
			if snp.POS >= regs["start"][i] and snp.POS<=regs["end"][i]:
				is_in_region=1
		if is_in_region==1:
			good_snps.append(snp)
			is_in_region=0
	return good_snps


def common_2snps(a,b):
	return [x for x in a if x in b]
def common_3snps(a,b,c):
	return [x for x in a if x in b and x in c]
def info_cover_filter(sam, thresh):
	vcf_filtered=[x for x in sam if x.INFO["DP"]>=thresh]
	return vcf_filtered

def read_vcf(name):
    vcf_raw=[]
    vcf_reader= vcf.Reader(open(name,"r"))
    for rec in vcf_reader:
        vcf_raw.append(rec)
    return vcf_raw
def write_vcf(name1,name2,good_vcf):
    writer=vcf.Writer(open(name1,"w"), vcf.Reader(open(name2,"r")))
    for i in range(0, len(good_vcf)):
        writer.write_record(good_vcf[i])
    writer.close()

def cover_filter(vcf_cont,cov_thresh):
    good_vcf=[]
    for rec in vcf_cont:
        i=0
        for sample in rec.samples: 
            if sample["DP"] < cov_thresh:
                i+=1
        if i==0:
            good_vcf.append(rec)
    return good_vcf
def qual_filter(vcf_cont, qual_thresh):
    good_vcf=[]
    for rec in vcf_cont:
        if rec.QUAL>=qual_thresh:
            good_vcf.append(rec)
    return good_vcf

