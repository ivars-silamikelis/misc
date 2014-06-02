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

def het_filter(vcf_raw,het_thresh):
    
    vcf_cont=[]    
    for rec in vcf_raw:
        i=0
        j=0
        for sample in rec.samples:
            if sample.data[0]!=None:
                if sample.data[0]=="1/1" or sample.data[0]=="0/1" or sample.data[0]=="0/0":
                    het=float(sample.data[1][0])/(sample.data[1][0]+sample.data[1][1])
                    if het<0.5:
                        het=1-het
                elif sample.data[0]=="2/2" or sample.data[0]=="0/2":
                    het=float(sample.data[1][0])/(sample.data[1][0]+sample.data[1][2])
                    if het<0.5:
                        het=1-het
                elif sample.data[0]=="1/2":
                    het=float(sample.data[1][1])/(sample.data[1][2]+sample.data[1][1])
                    if het<0.5:
                        het=1-het
                elif sample.data[0]=="1/3":
                    het=float(sample.data[1][1])/(sample.data[1][3]+sample.data[1][1])
                    if het<0.5:
                        het=1-het
                elif sample.data[0]=="2/3":
                    het=float(sample.data[1][2])/(sample.data[1][3]+sample.data[1][2])
                    if het<0.5:
                        het=1-het
                elif sample.data[0]=="3/3" or sample.data[0]=="0/3":
                    het=float(sample.data[1][0])/(sample.data[1][0]+sample.data[1][3])
                    if het<0.5:
                        het=1-het
                if het<het_thresh:
                    j+=1
            else:
                j+=1

        if i==0 and j==0:
            vcf_cont.append(rec)
    return vcf_cont


#visas pozicijas ar klaster snpiem. Ja kadam paraugam snp ir klasteri, tad abas pozicijas tiek ieliktas saraksta

def clst_filter(vcf_cont,clst_thresh): 
    pos=dict()
    for sample in vcf_cont[0].samples:
        pos[sample.sample]=[]	#dictionary inicializacija
    for rec in vcf_cont:
        for sample in rec.samples:
            if sample.data[0]!=None:            
                if sample.data[0]!="0/0":
                    pos[sample.sample].append(rec.POS)
    tmp=[]
    for sample in pos.keys():
        if len(pos[sample])>0:
            if np.abs(pos[sample][0]-pos[sample][1])<=clst_thresh:	#parbauda pirmo poziciju
               tmp.append(pos[sample][0])
            for i in range(1,len(pos[sample])-1):
                if np.abs(pos[sample][i-1]-pos[sample][i])<=clst_thresh or np.abs(pos[sample][i]-pos[sample][i+1])<=clst_thresh:
                    tmp.append(pos[sample][i])
            if np.abs(pos[sample][-1]-pos[sample][-2])<=clst_thresh:	#parbauda pedejo poziciju
                tmp.append(pos[sample][-1])
    bad_pos=sorted(set(tmp))
    good_vcf=[]
    i=0
    for rec in vcf_cont:
        if not rec.POS in bad_pos:
            good_vcf.append(rec)
    return good_vcf

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

