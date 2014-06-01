import vcf
import numpy as np

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
    i=0
    good_vcf=[]
    for rec in vcf_cont:
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
#writer=vcf.Writer(open("act.raw.hom"+het_thresh+".near"+clst_thresh+".vcf","w"), vcf.Reader(open("act.raw.vcf","r")))
#for i in range(0, len(good_vcf)):
#    writer.write_record(good_vcf[i])
#writer.close()
