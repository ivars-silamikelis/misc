import vcf
import filter_vcf as flt
import plot_vcf as plt_vcf
name="mtb_h37rv.raw.hom90.near10"
vcf_reader= vcf.Reader(open(name+".vcf","r"))
vcf_cont=[]
for rec in vcf_reader:
    vcf_cont.append(rec)

vcf_cov10=flt.cover_filter(vcf_cont, 10)
writer=vcf.Writer(open(name+".cov10.vcf","w"), vcf.Reader(open(name+".vcf","r")))
for i in range(0, len(vcf_cov10)):
    writer.write_record(vcf_cov10[i])
writer.close()
