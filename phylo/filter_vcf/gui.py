import Tkinter, tkFileDialog, Tkconstants

import filter_vcf as flt
import vcf
import numpy as np
vcf_raw=[]
filename=[]
top=Tkinter.Tk()
def callback():
    tkFileDialog.askopenfile(mode="r")
    filename.append(tkFileDialog.askopenfilename())
    print filename[0]+" loaded"
    vcf_reader= vcf.Reader(open(filename[0],"r"))
    for rec in vcf_reader:
        vcf_raw.append(rec)

def run_filter(vcf_r, t1, t2, t3, t4):
    
    t1=float(t1)
    t2=int(t2)
    t3=int(t3)
    t4=int(t4)
    good_vcf=flt.clst_filter(flt.qual_filter(flt.cover_filter(flt.het_filter(vcf_raw,t1),t3),t4),t2)
    #good_vcf_1=flt.het_filter(vcf_r,float(t1))
    #good_vcf_2=flt.clst_filter(good_vcf_1,float(t2))
    #good_vcf_3=flt.cover_filter(good_vcf_2,float(t3))
    #good_vcf=flt.qual_filter(good_vcf_3,float(t4))
    writer=vcf.Writer(open(filename[0]+".hom"+str(int(t1*100))+".near"+str(t2)+".cov"+str(t3)+".qual"+str(t4)+".vcf","w"), vcf.Reader(open(filename[0],"r")))
    for i in range(0, len(good_vcf)):
        writer.write_record(good_vcf[i])
    writer.close()
    print "File "+filename[0]+".hom"+str(int(t1*100))+".near"+str(t2)+".cov"+str(t3)+".qual"+str(t4)+".vcf created"
def gogogo():  
    run_filter(vcf_raw,het_flt_e.get(), clst_flt_e.get(), cov_flt_e.get(), qual_flt_e.get() )

bopen=Tkinter.Button(top,text="Open",command=callback)
bopen.pack()

Tkinter.Label(top, text="heterogenity filter").pack()
het_flt_e=Tkinter.Entry(top)
het_flt_e.pack()
het_flt_e.insert(0,0)

Tkinter.Label(top, text="close snp filter").pack()
clst_flt_e=Tkinter.Entry(top)
clst_flt_e.pack()
clst_flt_e.insert(0,0)

Tkinter.Label(top, text="coverage filter").pack()
cov_flt_e=Tkinter.Entry(top)
cov_flt_e.pack()
cov_flt_e.insert(0,0)


Tkinter.Label(top, text="quality filter").pack()
qual_flt_e=Tkinter.Entry(top)
qual_flt_e.pack()
qual_flt_e.insert(0,0)

brun=Tkinter.Button(top,text="Run",command=gogogo)
brun.pack()
top.mainloop
