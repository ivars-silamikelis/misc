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
	print "homozigosity filtering for" + name
	vs_homs=[x for x in varscan if x.INFO["HOM"]==1]
	#filtrs pec coverage un homozigotates
	
	gatk_homs=info_cover_filter(gatk,cover_thresh)
	sam_homs=info_cover_filter([x for x in sam if x.samples[0]["GT"]=="1/1"],cover_thresh)
	varscan_homs=cover_filter([x for x in varscan if x.INFO["HOM"]==1], cover_thresh)

	#starp snpcalleriem kopejie snp
	print "common snp filtering for" + name
	common_2gatk=common_2snps(gatk_homs,sam_homs)
	common_3gatk=common_3snps(gatk_homs,sam_homs,varscan_homs)
	#filtrs pec regioniem
	print "region filtering for" + name
	common_2gatk_final=keep_regs(filter_regs(common_2gatk, reg1),reg2)
	common_3gatk_final=keep_regs(filter_regs(common_3gatk, reg1),reg2)
	print "writing vcf for" + name
	write_vcf("filtered_gatk_sam/"+name+".gatk_sam.vcf",name+".gatk_ug.vcf",common_2gatk_final)
	write_vcf("filtered_gatk_sam_varscan/"+name+".gatk_sam_varscan.vcf",name+".gatk_ug.vcf",common_3gatk_final)

def filter_regs(snps, regs):
	#good_snps=snps
	bad_indexes=[]
	is_in_region=0	
	for snp_index in range(0, len(snps)): #good_snps:
		for i in range(0, len(regs["start"])):
			if snps[snp_index].POS >= regs["start"][i] and snps[snp_index].POS<=regs["end"][i]:
				bad_indexes.append(snp_index)
	for idx in bad_indexes[::-1]:
		print "removing "+ str(snps[idx].POS)
		snps.remove(snps[idx])
	return snps

def keep_regs(snps, regs):
	good_snps=[]
	is_in_region=0	
	for snp in snps:
		for i in range(0, len(regs["start"])):
			if snp.POS >= regs["start"][i] and snp.POS<=regs["end"][i]:
				is_in_region=1
		if is_in_region==1:
			print "keeping "+str(snp.POS)
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

